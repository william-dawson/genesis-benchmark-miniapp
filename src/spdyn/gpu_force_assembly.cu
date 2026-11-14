#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstddef>

extern "C" void genesis_nonbond_force_view(void **,void **);

namespace {
void checked(cudaError_t result) {
    if(result!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS force assembly: %s\n",cudaGetErrorString(result));
        std::abort();
    }
}
double *assembled=nullptr,*slabs=nullptr,*long_force=nullptr;
int *counts=nullptr,*starts=nullptr;
void *registered=nullptr;
double *mapped=nullptr;
std::size_t registered_bytes=0,output_bytes=0;
int saved_width=0,saved_cells=0,saved_capacity=0,saved_workers=0;
bool valid=false;
const double *bonded=nullptr;

__global__ void assemble(const double *cpu,const double *gpu,const double *fl,
                         const int *natom,const int *start,double *out,
                         int width,int cells,int capacity,int out_cells,int workers,int atoms,bool merge,const double *bond) {
    const int index=blockIdx.x*blockDim.x+threadIdx.x;
    const int cell=index/width,atom=index%width;
    if(cell>=out_cells)return;
    const std::size_t slot=3ull*(atom+width*cell);
    if(cell>=cells||atom>=natom[cell]) {
        out[slot]=0;out[slot+1]=0;out[slot+2]=0;return;
    }
    const std::size_t stride=3ull*width*capacity;
    double x=bond?bond[slot]:cpu[slot],y=bond?bond[slot+1]:cpu[slot+1],z=bond?bond[slot+2]:cpu[slot+2];
    if(merge) {x+=fl[slot];y+=fl[slot+1];z+=fl[slot+2];}
    const int dense=start[cell]+atom;
    x+=gpu[dense];y+=gpu[atoms+dense];z+=gpu[2*atoms+dense];
    for(int worker=1;!bond&&worker<workers;++worker) {
        x+=cpu[slot+stride*worker];y+=cpu[slot+1+stride*worker];z+=cpu[slot+2+stride*worker];
    }
    out[slot]=x;out[slot+1]=y;out[slot+2]=z;
}
}

extern "C" void genesis_force_assemble(const double *cpu,const double *fl,double *output,
    const int *natom,const int *start,int width,int cells,int capacity,int out_cells,int workers,int atoms,int merge,int mode,int download) {
    if(width<=0||cells<=0||capacity<cells||out_cells<cells||workers<=0||(mode!=1&&mode!=2))std::abort();
    const std::size_t slab_bytes=3ull*width*capacity*workers*sizeof(double);
    const std::size_t bytes=3ull*width*out_cells*sizeof(double);
    if(!assembled) {
        saved_width=width;saved_cells=cells;saved_capacity=capacity;saved_workers=workers;output_bytes=bytes;
        checked(cudaMalloc(&assembled,bytes));checked(cudaMalloc(&long_force,bytes));
        checked(cudaMalloc(&counts,cells*sizeof(int)));checked(cudaMalloc(&starts,cells*sizeof(int)));
        if(mode==1)checked(cudaMalloc(&slabs,slab_bytes));
    }
    if(saved_width!=width||saved_cells!=cells||saved_capacity!=capacity||saved_workers!=workers||output_bytes!=bytes) {
        std::fprintf(stderr,"GENESIS force assembly shape changed\n");std::abort();
    }
    void *raw_force=nullptr,*raw_stream=nullptr;
    genesis_nonbond_force_view(&raw_force,&raw_stream);
    if(!raw_force)std::abort();
    cudaStream_t stream=reinterpret_cast<cudaStream_t>(raw_stream);
    checked(cudaMemcpyAsync(counts,natom,cells*sizeof(int),cudaMemcpyHostToDevice,stream));
    checked(cudaMemcpyAsync(starts,start,cells*sizeof(int),cudaMemcpyHostToDevice,stream));
    const double *source=nullptr;
    if(bonded) { source=nullptr; } else if(mode==1) {
        checked(cudaMemcpyAsync(slabs,cpu,slab_bytes,cudaMemcpyHostToDevice,stream));source=slabs;
    } else {
        if(registered!=cpu||registered_bytes!=slab_bytes) {
            if(registered)checked(cudaHostUnregister(registered));
            registered=const_cast<double*>(cpu);registered_bytes=slab_bytes;
            checked(cudaHostRegister(registered,slab_bytes,cudaHostRegisterMapped));
            checked(cudaHostGetDevicePointer(&mapped,registered,0));
        }
        source=mapped;
    }
    if(merge)checked(cudaMemcpyAsync(long_force,fl,bytes,cudaMemcpyHostToDevice,stream));
    assemble<<<(width*out_cells+255)/256,256,0,stream>>>(source,static_cast<double*>(raw_force),long_force,
        counts,starts,assembled,width,cells,capacity,out_cells,workers,atoms,merge!=0,bonded);
    checked(cudaGetLastError());
    if(download)checked(cudaMemcpyAsync(output,assembled,bytes,cudaMemcpyDeviceToHost,stream));
    checked(cudaStreamSynchronize(stream));valid=true;
}

extern "C" int genesis_force_copy(double *destination,int width,int cells) {
    if(!valid)return 0;
    const std::size_t bytes=3ull*width*cells*sizeof(double);
    if(width!=saved_width||bytes>output_bytes)std::abort();
    checked(cudaMemcpy(destination,assembled,bytes,cudaMemcpyDeviceToDevice));
    return 1;
}

extern "C" void genesis_bonded_register(const double *values) { bonded=values; }
