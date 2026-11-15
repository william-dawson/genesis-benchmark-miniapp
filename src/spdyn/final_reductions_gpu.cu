#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
void checked(cudaError_t code) {
    if(code!=cudaSuccess){std::fprintf(stderr,"GENESIS final sums: %s\n",cudaGetErrorString(code));std::abort();}
}
__global__ void ordered_cells(const double *input,double *output,int cells,const int *failure) {
    __shared__ double tile[9*256];
    const int lane=threadIdx.x;
    double total=0;
    for(int start=0;start<cells;start+=256) {
        const int count=min(256,cells-start);
        for(int i=lane;i<9*count;i+=256)tile[i]=input[9*start+i];
        __syncthreads();
        // Each field keeps the original increasing-cell summation order.
        if(lane<9)for(int cell=0;cell<count;++cell)total=__dadd_rn(total,tile[9*cell+lane]);
        __syncthreads();
    }
    if(lane<9)output[lane]=total;
    if(!lane)output[9]=failure?*failure:0;
}
__device__ double warp_sum(double value) {
    for(int offset=16;offset;offset/=2)value=__dadd_rn(value,__shfl_down_sync(0xffffffff,value,offset));
    return value;
}
__global__ void tree_cells(const double *input,double *output,int cells,const int *failure) {
    __shared__ double sums[8];
    const int lane=threadIdx.x%32,warp=threadIdx.x/32,field=blockIdx.x;
    double value=0;
    for(int cell=threadIdx.x;cell<cells;cell+=256)value=__dadd_rn(value,input[9*cell+field]);
    value=warp_sum(value);
    if(!lane)sums[warp]=value;
    __syncthreads();
    if(!warp) {
        value=lane<8?sums[lane]:0;
        value=warp_sum(value);
        if(!lane){output[field]=value;if(!field)output[9]=failure?*failure:0;}
    }
}

}
extern "C" int genesis_final_reductions_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_FINAL_REDUCTIONS");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        if(end==text||*end||value<0||value>2){std::fprintf(stderr,"GENESIS_FINAL_REDUCTIONS must be 0, 1 or 2\n");std::abort();}
        return int(value);
    }();return mode;
}
extern "C" int genesis_final_reductions_tree() {
    static const int result=[] {
        const char *text=std::getenv("GENESIS_FINAL_REDUCTIONS_TREE");if(!text)return 0;
        if(!std::strcmp(text,"0"))return 0;
        if(!std::strcmp(text,"1"))return 1;
        std::fprintf(stderr,"GENESIS_FINAL_REDUCTIONS_TREE must be 0 or 1\n");std::abort();
    }();return result;
}
namespace {
void reduce(const double *input,double *output,int cells,const int *failure,int *bad) {
    if(cells<=0||!input||!output||!genesis_final_reductions_mode())std::abort();
    static double *device=nullptr,*host=nullptr;
    if(!device){checked(cudaMalloc(&device,10*sizeof(double)));checked(cudaMallocHost(&host,10*sizeof(double)));}
    if(genesis_final_reductions_tree())tree_cells<<<9,256>>>(input,device,cells,failure);
    else ordered_cells<<<1,256>>>(input,device,cells,failure);
    checked(cudaGetLastError());
    checked(cudaMemcpyAsync(host,device,(failure?10:9)*sizeof(double),cudaMemcpyDeviceToHost,0));
    checked(cudaStreamSynchronize(0));std::memcpy(output,host,9*sizeof(double));
    if(bad)*bad=static_cast<int>(host[9]);
}
}
extern "C" void genesis_final_reduce(const double *input,double *output,int cells) {
    reduce(input,output,cells,nullptr,nullptr);
}
extern "C" void genesis_final_reduce_checked(const double *input,double *output,int cells,const int *failure,int *bad) {
    if(!failure||!bad)std::abort();
    reduce(input,output,cells,failure,bad);
}
