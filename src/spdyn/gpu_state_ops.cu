#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>

namespace {
__global__ void remove_translation(double *velocity,const int *natom,
                                   int width,int cells,double x,double y,double z) {
    const int index=blockIdx.x*blockDim.x+threadIdx.x;
    const int cell=index/width,atom=index%width;
    if(cell>=cells||atom>=natom[cell])return;
    const size_t offset=3ull*index;
    velocity[offset]-=x;velocity[offset+1]-=y;velocity[offset+2]-=z;
}
}

extern "C" void genesis_stop_translation(double *velocity,const int *natom,
    int width,int cells,double x,double y,double z) {
    if(width<=0||cells<=0)std::abort();
    remove_translation<<<(width*cells+255)/256,256,0,0>>>(velocity,natom,width,cells,x,y,z);
    const cudaError_t status=cudaGetLastError();
    if(status!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS translation removal: %s\n",cudaGetErrorString(status));
        std::abort();
    }
}

extern "C" void genesis_poison_work(double *first,double *second,int width,int cells) {
    const size_t bytes=3ull*width*cells*sizeof(double);
    for(double *values:{first,second}) {
        const cudaError_t status=cudaMemsetAsync(values,0xff,bytes,0);
        if(status!=cudaSuccess) {
            std::fprintf(stderr,"GENESIS scratch poison: %s\n",cudaGetErrorString(status));
            std::abort();
        }
    }
}
