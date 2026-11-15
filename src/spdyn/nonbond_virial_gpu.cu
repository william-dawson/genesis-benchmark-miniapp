#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>

extern "C" int genesis_coord_active();
namespace {
bool requested=false,needed=false,ready=false,consumer_pending=false;
int atoms=0;
double *partials=nullptr,*totals=nullptr,*host=nullptr;
cudaEvent_t complete;
void require(bool ok,const char *message) {
    if(!ok) {std::fprintf(stderr,"GENESIS nonbond virial: %s\n",message);std::abort();}
}
void checked(cudaError_t code) {require(code==cudaSuccess,cudaGetErrorString(code));}
bool selected(const char *name,const char *expected) {
    const char *text=std::getenv(name);return text&&std::strcmp(text,expected)==0;
}
int mode() {
    static const int result=[] {
        const char *text=std::getenv("GENESIS_NONBOND_VIRIAL");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");return int(value);
    }();return result;
}
__global__ void dot_partials(const double *x,const double *force,double *output,int count) {
    __shared__ double sums[3][256];
    const int lane=threadIdx.x,index=blockIdx.x*256+lane;
    for(int axis=0;axis<3;++axis)
        sums[axis][lane]=index<count?__dmul_rn(x[axis*count+index],force[axis*count+index]):0.0;
    __syncthreads();
    for(int offset=128;offset;offset/=2) {
        if(lane<offset)for(int axis=0;axis<3;++axis)
            sums[axis][lane]=__dadd_rn(sums[axis][lane],sums[axis][lane+offset]);
        __syncthreads();
    }
    if(!lane)for(int axis=0;axis<3;++axis)output[3*blockIdx.x+axis]=sums[axis][0];
}
__global__ void dot_finish(const double *parts,double *output,int count) {
    __shared__ double sums[3][256];
    const int lane=threadIdx.x;
    for(int axis=0;axis<3;++axis) {
        double value=0;
        for(int index=lane;index<count;index+=256)value=__dadd_rn(value,parts[3*index+axis]);
        sums[axis][lane]=value;
    }
    __syncthreads();
    for(int offset=128;offset;offset/=2) {
        if(lane<offset)for(int axis=0;axis<3;++axis)
            sums[axis][lane]=__dadd_rn(sums[axis][lane],sums[axis][lane+offset]);
        __syncthreads();
    }
    if(!lane)for(int axis=0;axis<3;++axis)output[axis]=sums[axis][0];
}
}
extern "C" void genesis_nonbond_virial_request(int eligible,int count,int need) {
    require(!requested&&!ready&&!consumer_pending,"previous result was not consumed");
    requested=mode()&&eligible&&selected("GENESIS_FORCE_ASSEMBLY","2")&&selected("GENESIS_FORCE_BONDED","1")
        &&!selected("GENESIS_STEP_VERIFY","1")&&!selected("GENESIS_STEP_VERIFY","2");
    needed=need!=0;
    if(!requested)return;
    require(count>0,"invalid atom count");
    if(atoms)require(count==atoms,"atom count changed");
    atoms=count;
}
extern "C" int genesis_nonbond_virial_host_required() {
    return !(requested&&mode()==1&&genesis_coord_active()&&!selected("GENESIS_FORCE_VERIFY","1")
        &&!selected("GENESIS_COORD_VERIFY","1"));
}
extern "C" int genesis_nonbond_virial_stage(const void *coordinates,const void *force,size_t bytes,int width,void *queue,double *host_coordinates,double *host_force) {
    if(!requested)return 0;
    requested=false;
    require(width==sizeof(double)&&bytes==3ull*atoms*sizeof(double),"invalid force shape");
    if(!needed) {
        if(selected("GENESIS_FORCE_VERIFY","1"))return 0;
        consumer_pending=true;
        if(mode()==2) {
            for(int i=0;i<3*atoms;++i)host_coordinates[i]=host_force[i]=std::numeric_limits<double>::quiet_NaN();
            std::printf("GENESIS_NONBOND_VIRIAL: no-virial mirrors poisoned\n");
        }
        return 1;
    }
    const int blocks=(atoms+255)/256;
    if(!partials) {
        checked(cudaMalloc(&partials,3ull*blocks*sizeof(double)));checked(cudaMalloc(&totals,3*sizeof(double)));
        checked(cudaMallocHost(&host,3*sizeof(double)));checked(cudaEventCreateWithFlags(&complete,cudaEventDisableTiming));
    }
    const auto stream=reinterpret_cast<cudaStream_t>(queue);
    dot_partials<<<blocks,256,0,stream>>>(static_cast<const double*>(coordinates),static_cast<const double*>(force),partials,atoms);
    checked(cudaGetLastError());
    dot_finish<<<1,256,0,stream>>>(partials,totals,blocks);checked(cudaGetLastError());
    checked(cudaMemcpyAsync(host,totals,3*sizeof(double),cudaMemcpyDeviceToHost,stream));
    checked(cudaEventRecord(complete,stream));ready=true;
    const bool retire=mode()==1&&!selected("GENESIS_FORCE_VERIFY","1");
    if(retire)consumer_pending=true;
    return retire;
}
extern "C" int genesis_nonbond_virial_get(double *output) {
    if(!ready)return 0;
    checked(cudaEventSynchronize(complete));ready=false;
    for(int axis=0;axis<3;++axis) {require(std::isfinite(host[axis]),"nonfinite result");output[axis]=host[axis];}
    return mode();
}
extern "C" void genesis_nonbond_virial_poison(double *coordinates,double *force) {
    if(selected("GENESIS_FORCE_VERIFY","1"))return;
    consumer_pending=true;
    for(int i=0;i<3*atoms;++i)coordinates[i]=force[i]=std::numeric_limits<double>::quiet_NaN();
    std::printf("GENESIS_NONBOND_VIRIAL: retired coordinate and force mirrors poisoned\n");
}

extern "C" void genesis_nonbond_virial_consume() {
    if(!consumer_pending)return;
    consumer_pending=false;
    if(mode()==2)std::printf("GENESIS_NONBOND_VIRIAL: direct GPU force consumer pass\n");
}
