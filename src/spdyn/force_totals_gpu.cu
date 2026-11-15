#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
void checked(cudaError_t status) {
    if(status!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS bonded totals: %s\n",cudaGetErrorString(status));
        std::abort();
    }
}

__device__ double add_cells(double total,const double *parts,int field,
                           int worker,int cells,int workers) {
    for(int cell=worker;cell<cells;cell+=workers)
        total=__dadd_rn(total,parts[field+10*cell]);
    return total;
}

__global__ void finish_totals(const double *bond,const double *angle,
    const double *dihedral,const double *improper,const double *nb14,
    const double *correction,double *output,int cells,int workers) {
    const int index=blockIdx.x*blockDim.x+threadIdx.x;
    if(index>=10*workers) return;
    const int field=index%10,worker=index/10;
    // Preserve the host's term/cell order within each worker accumulator.
    double total=add_cells(0.0,bond,field,worker,cells,workers);
    total=add_cells(total,angle,field,worker,cells,workers);
    total=add_cells(total,dihedral,field,worker,cells,workers);
    total=add_cells(total,improper,field,worker,cells,workers);
    total=add_cells(total,nb14,field,worker,cells,workers);
    total=add_cells(total,correction,field,worker,cells,workers);
    output[index]=total;
}

__global__ void staged_totals(const double *bond,const double *angle,
    const double *dihedral,const double *improper,const double *nb14,
    const double *correction,double *output,int cells,int workers,int slots) {
    extern __shared__ double parts[];
    const int worker=blockIdx.x,plane=10*slots;
    for(int index=threadIdx.x;index<plane;index+=blockDim.x) {
        const int cell=worker+(index/10)*workers,field=index%10;
        if(cell<cells) {
            const int source=field+10*cell;
            parts[index]=bond[source];
            parts[plane+index]=angle[source];
            parts[2*plane+index]=dihedral[source];
            parts[3*plane+index]=improper[source];
            parts[4*plane+index]=nb14[source];
            parts[5*plane+index]=correction[source];
        }
    }
    __syncthreads();
    if(threadIdx.x<10) {
        const int field=threadIdx.x,count=(cells-worker+workers-1)/workers;
        double total=0.0;
        // Cooperative loads hide latency; shared-memory sums retain CPU order.
        for(int term=0;term<6;++term)
            for(int cell=0;cell<count;++cell)
                total=__dadd_rn(total,parts[term*plane+10*cell+field]);
        output[field+10*worker]=total;
    }
}
}

extern "C" int genesis_bonded_totals_mode() {
    static const int mode=[]() {
        const char *text=std::getenv("GENESIS_BONDED_TOTALS");
        if(!text) return 0;
        char *end=nullptr;
        const long value=std::strtol(text,&end,10);
        if(end==text||*end||value<0||value>2) {
            std::fprintf(stderr,"GENESIS_BONDED_TOTALS must be 0, 1 or 2\n");
            std::abort();
        }
        return int(value);
    }();
    return mode;
}

extern "C" void genesis_bonded_totals(const double *bond,const double *angle,
    const double *dihedral,const double *improper,const double *nb14,
    const double *correction,int cells,int workers,cudaStream_t stream,double *output) {
    if(cells<=0||workers<=0) std::abort();
    static double *device=nullptr,*host=nullptr;
    static int capacity=0;
    const int count=10*workers;
    const size_t bytes=size_t(count)*sizeof(double);
    if(count>capacity) {
        if(device) {checked(cudaFree(device));checked(cudaFreeHost(host));}
        checked(cudaMalloc(&device,bytes));
        checked(cudaMallocHost(&host,bytes));
        capacity=count;
    }
    static const int staged=[]() {
        const char *text=std::getenv("GENESIS_BONDED_TOTALS_STAGED");
        if(!text) return 1;
        if(!std::strcmp(text,"0")) return 0;
        if(!std::strcmp(text,"1")) return 1;
        std::fprintf(stderr,"GENESIS_BONDED_TOTALS_STAGED must be 0 or 1\n");
        std::abort();
    }();
    const int slots=(cells+workers-1)/workers;
    const size_t shared_bytes=60ull*slots*sizeof(double);
    if(staged&&shared_bytes<=48*1024)
        staged_totals<<<workers,128,shared_bytes,stream>>>(bond,angle,dihedral,
            improper,nb14,correction,device,cells,workers,slots);
    else
        finish_totals<<<(count+127)/128,128,0,stream>>>(bond,angle,dihedral,
            improper,nb14,correction,device,cells,workers);
    checked(cudaGetLastError());
    checked(cudaMemcpyAsync(host,device,bytes,cudaMemcpyDeviceToHost,stream));
    checked(cudaStreamSynchronize(stream));
    std::memcpy(output,host,bytes);
}
