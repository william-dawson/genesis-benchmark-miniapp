#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include "../config.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <chrono>

extern "C" void genesis_hgr_scatter(double *,double *,int *,int *,int *,int *,
    const double *,const int *,const int *,int,int,int,int,int);
namespace {
void require(bool value,const char *message) {
    if(!value) {std::fprintf(stderr,"HGr migration: %s\n",message);std::abort();}
}
void checked(cudaError_t code) {
    if(code!=cudaSuccess) {std::fprintf(stderr,"HGr CUDA: %s\n",cudaGetErrorString(code));std::abort();}
}
int flag(const char *name) {
    const char *value=std::getenv(name);
    if(!value)return 0;
    require(!std::strcmp(value,"0")||!std::strcmp(value,"1"),name);
    return value[0]=='1';
}
struct Attributes {double charge,mass;int cls,solute;};
struct State {
    const double *x=nullptr,*v=nullptr,*dist=nullptr;
    const int *groups=nullptr,*counts=nullptr;
    int width=0,cells=0,group_width=0,sites=0,connect=0,capacity=0,globals=0,total=0,atoms=0;
    int *keys[2]{},*values[2]{},*prefix=nullptr,*ids=nullptr,*ranges=nullptr,*error=nullptr,*integers=nullptr;
    std::int16_t *map=nullptr;
    double *real=nullptr,*host_real=nullptr;
    int *host_integer=nullptr,*host_ranges=nullptr;
    Attributes *attributes=nullptr;
    void *temporary=nullptr;size_t temporary_bytes=0;
    cudaEvent_t ready=nullptr,start=nullptr,finish=nullptr;
    bool latched=false,timing=false;
} state;
__device__ int component(double r,double origin,double box,double cell,int count) {
    const double shifted=__dadd_rn(r,-origin);
    const double offset=__dadd_rn(__dmul_rn(box,0.5),-__dmul_rn(box,round(shifted/box)));
    int result=int(__dadd_rn(shifted,offset)/cell);
    if(result==count)--result;
    return result;
}
__global__ void hgr_classify(const double *x,const int *groups,const int *counts,const int *prefix,
    const std::int16_t *map,int *keys,int *values,int *error,int width,int cells,int group_width,int sites,int connect,
    double3 box,double3 origin,double3 size,int3 grid) {
    const int slot=blockIdx.x*blockDim.x+threadIdx.x;
    const int bucket=slot/group_width,entry=slot%group_width,cell=bucket/connect,kind=bucket%connect;
    if(cell>=cells||entry>=counts[bucket])return;
    const int output=prefix[bucket]+entry;
    keys[output]=0;values[output]=slot;
    const int atom=groups[size_t(sites)*slot]-1;
    if(atom<0||atom>=width) {atomicExch(error,1);return;}
    const double *r=x+3ull*(atom+width*cell);
    const int cx=component(r[0],origin.x,box.x,size.x,grid.x);
    const int cy=component(r[1],origin.y,box.y,size.y,grid.y);
    const int cz=component(r[2],origin.z,box.z,size.z,grid.z);
    if(cx<0||cx>=grid.x||cy<0||cy>=grid.y||cz<0||cz>=grid.z) {atomicExch(error,2);return;}
    const int destination=int(map[cx+grid.x*(cy+grid.y*cz)])-1;
    if(destination<0||destination>=cells) {atomicExch(error,3);return;}
    keys[output]=kind+connect*(destination==cell?cell:cells+destination);
}
__global__ void hgr_pack(const double *x,const double *v,const double *dist,const int *groups,const int *ids,
    const Attributes *attributes,const int *keys,const int *values,int *ranges,double *real,int *integer,int *error,
    int width,int group_width,int sites,int connect,int total,int atoms) {
    const int item=blockIdx.x*blockDim.x+threadIdx.x;
    if(item>=total)return;
    const int key=keys[item],slot=values[item],bucket=slot/group_width;
    const int cell=bucket/connect,kind=bucket%connect;
    if(item==0||keys[item-1]!=key)ranges[2*key]=item;
    if(item==total-1||keys[item+1]!=key)ranges[2*key+1]=item+1;
    for(int list=0;list<sites;++list) {
        double *r=real+9ull*(sites*item+list);
        int *i=integer+3ull*(sites*item+list);
        for(int c=0;c<9;++c)r[c]=0;
        i[0]=i[1]=i[2]=0;
        if(list>kind+1)continue;
        const size_t group=size_t(sites)*slot+list;
        const int atom=groups[group]-1;
        if(atom<0||atom>=width) {atomicExch(error,4);continue;}
        const size_t input=atom+size_t(width)*cell;
        const int id=ids[input];
        if(id<=0||id>atoms) {atomicExch(error,5);continue;}
        const Attributes a=attributes[id];
        for(int c=0;c<3;++c) {r[c]=x[3*input+c];r[3+c]=v[3*input+c];}
        r[6]=a.charge;r[7]=a.mass;r[8]=dist[group];
        i[0]=a.cls;i[1]=id;i[2]=a.solute;
    }
}
}
extern "C" int genesis_hgr_enabled() {
    static int value=flag("GENESIS_HGR_GPU");
#if defined(_SINGLE) || defined(_MIXED)
    require(!value,"requires double precision");
#endif
    return value;
}
extern "C" int genesis_hgr_verify() {static int value=flag("GENESIS_HGR_VERIFY");return value;}
extern "C" void genesis_hgr_latch(const double *x,const double *v,const int *groups,const int *counts,
    const double *dist,int width,int cells,int group_width,int sites,int connect) {
    if(!genesis_hgr_enabled())return;
    auto &s=state;
    if(!s.ready)checked(cudaEventCreateWithFlags(&s.ready,cudaEventDisableTiming));
    if(s.width)require(s.width==width&&s.cells==cells&&s.group_width==group_width&&s.sites==sites&&s.connect==connect,"shape changed");
    s.x=x;s.v=v;s.groups=groups;s.counts=counts;s.dist=dist;
    s.width=width;s.cells=cells;s.group_width=group_width;s.sites=sites;s.connect=connect;
    checked(cudaEventRecord(s.ready,0));s.latched=true;
}
extern "C" int genesis_hgr_outgoing(const int *counts,const int *groups,const int *ids,const int *solute,
    const int *cls,const double *charge,const double *mass,const std::int16_t *map,
    const double *box,const double *origin,const double *size,const int *grid,int atoms,
    int *stay,int *move,double *stay_real,double *move_real,int *stay_integer,int *move_integer,
    int stay_capacity,int move_capacity) {
    auto &s=state;
    if(!genesis_hgr_enabled()||!s.latched)return 0;
    s.latched=false;
    auto begin=std::chrono::steady_clock::now();
    const int buckets=s.connect*s.cells,globals=grid[0]*grid[1]*grid[2];
    std::vector<int> prefix(buckets+1);
    for(int b=0;b<buckets;++b) {
        require(counts[b]>=0&&counts[b]<=s.group_width,"invalid group count");
        prefix[b+1]=prefix[b]+counts[b];
    }
    s.total=prefix.back();require(s.total>0&&atoms>0,"empty groups");
    if(!s.capacity) {
        s.capacity=s.total;s.globals=globals;s.atoms=atoms;
        for(auto &p:s.keys)checked(cudaMalloc(&p,size_t(s.capacity)*sizeof(int)));
        for(auto &p:s.values)checked(cudaMalloc(&p,size_t(s.capacity)*sizeof(int)));
        checked(cudaMalloc(&s.prefix,size_t(buckets+1)*sizeof(int)));
        checked(cudaMalloc(&s.ids,size_t(s.width)*s.cells*sizeof(int)));
        checked(cudaMalloc(&s.map,size_t(globals)*sizeof(std::int16_t)));
        checked(cudaMalloc(&s.ranges,4ull*buckets*sizeof(int)));
        checked(cudaMalloc(&s.error,sizeof(int)));
        checked(cudaMalloc(&s.real,9ull*s.sites*s.capacity*sizeof(double)));
        checked(cudaMalloc(&s.integers,3ull*s.sites*s.capacity*sizeof(int)));
        checked(cudaMallocHost(&s.host_real,9ull*s.sites*s.capacity*sizeof(double)));
        checked(cudaMallocHost(&s.host_integer,3ull*s.sites*s.capacity*sizeof(int)));
        checked(cudaMallocHost(&s.host_ranges,4ull*buckets*sizeof(int)));
        checked(cub::DeviceRadixSort::SortPairs(nullptr,s.temporary_bytes,s.keys[0],s.keys[1],s.values[0],s.values[1],s.capacity));
        checked(cudaMalloc(&s.temporary,s.temporary_bytes));
        std::vector<Attributes> attributes(size_t(atoms)+1);
        std::vector<bool> seen(size_t(atoms)+1);
        for(int b=0;b<buckets;++b)for(int entry=0;entry<counts[b];++entry)for(int l=0;l<b%s.connect+2;++l) {
            const int atom=groups[s.sites*(entry+s.group_width*b)+l]-1;
            require(atom>=0&&atom<s.width,"invalid initial atom");
            const size_t input=atom+size_t(s.width)*(b/s.connect);
            const int id=ids[input];require(id>0&&id<=atoms&&!seen[id],"duplicate or invalid initial global ID");
            seen[id]=true;attributes[id]={charge[input],mass[input],cls[input],solute[input]};
        }
        checked(cudaMalloc(&s.attributes,attributes.size()*sizeof(Attributes)));
        checked(cudaMemcpy(s.attributes,attributes.data(),attributes.size()*sizeof(Attributes),cudaMemcpyHostToDevice));
        s.timing=flag("GENESIS_HGR_TIMING");
        if(s.timing) {checked(cudaEventCreate(&s.start));checked(cudaEventCreate(&s.finish));}
    }
    require(globals==s.globals&&atoms==s.atoms&&s.total==s.capacity,"global shape or group count changed");
    checked(cudaStreamWaitEvent(0,s.ready,0));
    checked(cudaMemcpy(s.ids,ids,size_t(s.width)*s.cells*sizeof(int),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.prefix,prefix.data(),size_t(buckets+1)*sizeof(int),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.map,map,size_t(globals)*sizeof(std::int16_t),cudaMemcpyHostToDevice));
    checked(cudaMemset(s.ranges,0xff,4ull*buckets*sizeof(int)));checked(cudaMemset(s.error,0,sizeof(int)));
    if(s.timing)checked(cudaEventRecord(s.start,0));
    hgr_classify<<<(s.group_width*buckets+255)/256,256>>>(s.x,s.groups,s.counts,s.prefix,s.map,s.keys[0],s.values[0],s.error,
        s.width,s.cells,s.group_width,s.sites,s.connect,make_double3(box[0],box[1],box[2]),
        make_double3(origin[0],origin[1],origin[2]),make_double3(size[0],size[1],size[2]),make_int3(grid[0],grid[1],grid[2]));
    checked(cudaGetLastError());
    int bits=0;while((1ull<<bits)<2ull*buckets)++bits;
    checked(cub::DeviceRadixSort::SortPairs(s.temporary,s.temporary_bytes,s.keys[0],s.keys[1],s.values[0],s.values[1],s.total,0,bits));
    hgr_pack<<<(s.total+255)/256,256>>>(s.x,s.v,s.dist,s.groups,s.ids,s.attributes,s.keys[1],s.values[1],s.ranges,
        s.real,s.integers,s.error,s.width,s.group_width,s.sites,s.connect,s.total,s.atoms);
    checked(cudaGetLastError());if(s.timing)checked(cudaEventRecord(s.finish,0));
    checked(cudaMemcpy(s.host_real,s.real,9ull*s.sites*s.total*sizeof(double),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(s.host_integer,s.integers,3ull*s.sites*s.total*sizeof(int),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(s.host_ranges,s.ranges,4ull*buckets*sizeof(int),cudaMemcpyDeviceToHost));
    int error=0;checked(cudaMemcpy(&error,s.error,sizeof(int),cudaMemcpyDeviceToHost));require(!error,"invalid classification/atom");
    auto copied=std::chrono::steady_clock::now();
    for(int b=0;b<2*buckets;++b) {
        const int first=s.host_ranges[2*b],last=s.host_ranges[2*b+1],count=first<0?0:last-first;
        require(count>=0&&count<=(b<buckets?stay_capacity:move_capacity),"destination overflow");
    }
    genesis_hgr_scatter(stay_real,move_real,stay_integer,move_integer,stay,move,s.host_real,s.host_integer,
        s.host_ranges,s.cells,s.connect,s.sites,stay_capacity,move_capacity);
    if(s.timing) {
        float ms=0;checked(cudaEventElapsedTime(&ms,s.start,s.finish));
        std::printf("GENESIS_HGR_TIMING: %d %.9f %.9f %.9f\n",s.total,double(ms)*1e-3,
            std::chrono::duration<double>(copied-begin).count(),std::chrono::duration<double>(std::chrono::steady_clock::now()-copied).count());
    }
    return 1;
}
extern "C" void genesis_hgr_check(const int *stay,const int *move,const double *stay_real,const double *move_real,
    const int *stay_integer,const int *move_integer,int stay_capacity,int move_capacity) {
    const auto &s=state;const int buckets=s.cells*s.connect;int seen=0;
    for(int b=0;b<2*buckets;++b) {
        const int first=s.host_ranges[2*b],last=s.host_ranges[2*b+1],count=first<0?0:last-first;
        const int source=b%buckets,sites=source%s.connect+2,capacity=b<buckets?stay_capacity:move_capacity;
        require((b<buckets?stay:move)[source]==count,"CPU/GPU counts differ");
        for(int k=0;k<count;++k) {
            const size_t output=k+size_t(capacity)*source,input=first+k;
            require(!std::memcmp((b<buckets?stay_real:move_real)+9ull*s.sites*output,s.host_real+9ull*s.sites*input,
                9ull*sites*sizeof(double)),"CPU/GPU real records differ");
            require(!std::memcmp((b<buckets?stay_integer:move_integer)+3ull*s.sites*output,s.host_integer+3ull*s.sites*input,
                3ull*sites*sizeof(int)),"CPU/GPU identities/order differ");
        }
        seen+=count;
    }
    require(seen==s.total,"group conservation failed");
    std::printf("GENESIS_HGR_VERIFY: exact packing pass %d\n",seen);
}
