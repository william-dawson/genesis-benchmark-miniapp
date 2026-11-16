#include <cuda_runtime.h>
#include "precision_types.h"
using genesis::IntegrationReal;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace {
void checked(cudaError_t status) {
    if(status!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS long force: %s\n",cudaGetErrorString(status));
        std::abort();
    }
}
int mode() {
    static const int selected=[]() {
        const char *text=std::getenv("GENESIS_LONG_FORCE");
        if(!text) return 0;
        if(!std::strcmp(text,"0")) return 0;
        if(!std::strcmp(text,"1")) return 1;
        if(!std::strcmp(text,"2")) return 2;
        std::fprintf(stderr,"GENESIS_LONG_FORCE must be 0, 1 or 2\n");
        std::abort();
    }();
    return selected;
}
bool requested=false,valid=false;
IntegrationReal *padded=nullptr,*observed=nullptr,*integrator=nullptr,*published=nullptr;
int *counts=nullptr,*offsets=nullptr;
int saved_width=0,saved_cells=0;
std::vector<int> host_counts;

__global__ void scatter_long(const IntegrationReal *packed,IntegrationReal *output,
    const int *natom,const int *starts,int width,int cells) {
    const int index=blockIdx.x*blockDim.x+threadIdx.x;
    const int cell=index/width,atom=index%width;
    if(cell>=cells) return;
    const size_t destination=3ull*index;
    if(atom<natom[cell]) {
        const size_t source=3ull*(starts[cell]+atom);
        output[destination]=packed[source];
        output[destination+1]=packed[source+1];
        output[destination+2]=packed[source+2];
    } else {
        output[destination]=0;output[destination+1]=0;output[destination+2]=0;
    }
}
}

extern "C" void genesis_long_force_request(int eligible) {
    requested=false;
    valid=false;
    if(!eligible||mode()==0) return;
    static const bool compatible=[]() {
        const char *diagnostic=std::getenv("GENESIS_STEP_VERIFY");
        const bool reference=diagnostic&&(!std::strcmp(diagnostic,"1")||!std::strcmp(diagnostic,"2"));
        const char *pme=std::getenv("GENESIS_PME_GPU"),*packed=std::getenv("GENESIS_PME_PACKED");
        return !reference&&pme&&packed&&!std::strcmp(pme,"1")&&!std::strcmp(packed,"1");
    }();
    requested=compatible;
}
extern "C" int genesis_long_force_requested_() {return requested;}
extern "C" int genesis_long_force_reference() {return mode()==2;}
extern "C" int genesis_long_force_available() {return valid;}

extern "C" void genesis_long_force_publish(const IntegrationReal *force,const int *natom,
    const int *starts,int width,int cells,cudaStream_t stream) {
    if(!requested||width<=0||cells<=0) std::abort();
    if(!padded) {
        saved_width=width;saved_cells=cells;
        checked(cudaMalloc(&padded,3ull*width*cells*sizeof(IntegrationReal)));
        checked(cudaMalloc(&counts,cells*sizeof(int)));
        checked(cudaMalloc(&offsets,cells*sizeof(int)));
        if(mode()==2) checked(cudaMallocHost(&observed,3ull*width*cells*sizeof(IntegrationReal)));
    }
    if(width!=saved_width||cells!=saved_cells) std::abort();
    checked(cudaMemcpyAsync(counts,natom,cells*sizeof(int),cudaMemcpyHostToDevice,stream));
    checked(cudaMemcpyAsync(offsets,starts,cells*sizeof(int),cudaMemcpyHostToDevice,stream));
    static const bool direct=[]() {
        const char *text=std::getenv("GENESIS_LONG_DIRECT");
        if(!text||!std::strcmp(text,"1")) return true;
        if(!std::strcmp(text,"0")) return false;
        std::fprintf(stderr,"GENESIS_LONG_DIRECT must be 0 or 1\n");std::abort();
    }();
    published=direct&&integrator?integrator:padded;
    scatter_long<<<(width*cells+255)/256,256,0,stream>>>(force,published,counts,offsets,width,cells);
    checked(cudaGetLastError());
    if(mode()==2) host_counts.assign(natom,natom+cells);
    valid=true;
}

extern "C" int genesis_long_force_copy(IntegrationReal *destination,IntegrationReal *reference,int width,int cells) {
    if(!valid) return 0;
    if(width!=saved_width||cells!=saved_cells) std::abort();
    const size_t bytes=3ull*width*cells*sizeof(IntegrationReal);
    if(integrator&&integrator!=destination) std::abort();
    integrator=destination;
    // PME synchronizes publication; only the first handoff needs a copy.
    if(published!=destination) {
        checked(cudaMemcpy(destination,published,bytes,cudaMemcpyDeviceToDevice));
        checked(cudaStreamSynchronize(nullptr));
    }
    if(mode()==2) {
        checked(cudaMemcpy(observed,published,bytes,cudaMemcpyDeviceToHost));
        for(int cell=0;cell<cells;++cell) {
            const size_t first=3ull*width*cell,count=3ull*host_counts[cell];
            if(std::memcmp(observed+first,reference+first,count*sizeof(IntegrationReal))) {
                std::fprintf(stderr,"GENESIS long-force publication differs from CPU scatter\n");
                std::abort();
            }
            for(size_t atom=first;atom<first+count;++atom)
                reference[atom]=std::numeric_limits<IntegrationReal>::quiet_NaN();
        }
        std::puts("GENESIS_LONG_FORCE_VERIFY: exact publication and CPU mirror poisoned");
    }
    valid=false;
    return 1;
}
