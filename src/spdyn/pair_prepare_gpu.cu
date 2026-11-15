#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "migration_views.h"
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"pair preparation: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_incoming_pair_view(ParticleView*,const double**);
extern "C" int genesis_exclusion_prefix(const int**,const int**);
namespace {
ParticleView state{};const double *translation=nullptr,*host_coords=nullptr;const int *host_classes=nullptr,*counts=nullptr,*starts=nullptr;
bool active=false,prepared=false,fresh=false;const double *destination=nullptr;const int *class_destination=nullptr;cudaEvent_t ready=nullptr;int total=0;
__global__ void pack(ParticleView v,const double *translation,const int *counts,const int *starts,double *out,int *classes,int atoms){
 int c=blockIdx.x;for(int a=threadIdx.x;a<counts[c];a+=blockDim.x){int i=a+v.width*c,j=a+starts[c];
 for(int k=0;k<3;++k)out[j+size_t(k)*atoms]=__dadd_rn(v.x[3*i+k],translation[3*i+k]);out[j+size_t(3)*atoms]=v.q[i];classes[j]=v.classes[i];}
}
}
extern "C" int genesis_pair_prepare_mode(){static int mode=[](){const char *s=std::getenv("GENESIS_PAIR_PREPARE_GPU");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;if(!std::strcmp(s,"2"))return 2;std::abort();}();return mode;}
extern "C" int genesis_pair_prepare_begin(const double *coords,const int *classes,int cells,int width,int atoms){
 active=false;prepared=false;fresh=false;if(!genesis_pair_prepare_mode())return 0;
 if(!genesis_incoming_pair_view(&state,&translation) || !genesis_exclusion_prefix(&counts,&starts))return 0;
 if(state.cells!=cells || state.width!=width || atoms<=0)std::abort();host_coords=coords;host_classes=classes;total=atoms;active=true;return 1;
}
extern "C" int genesis_pair_prepare_copy(double *out,int *classes,int *out_counts,int *out_starts,const double *coords,const int *host_cls,int atoms,int ncell,cudaStream_t queue){
 if(!active)return 0;if(coords!=host_coords || host_cls!=host_classes || atoms!=total || ncell<state.cells){std::fprintf(stderr,"Pair preparation shape/identity mismatch: coordinates=%d classes=%d atoms=%d/%d cells=%d/%d\n",coords==host_coords,host_cls==host_classes,atoms,total,ncell,state.cells);std::abort();}
 pack<<<state.cells,128,0,queue>>>(state,translation,counts,starts,out,classes,atoms);CHECK(cudaGetLastError());
 CHECK(cudaMemcpyAsync(out_counts,counts,state.cells*4,cudaMemcpyDeviceToDevice,queue));CHECK(cudaMemcpyAsync(out_starts,starts,state.cells*4,cudaMemcpyDeviceToDevice,queue));
 if(ncell>state.cells){CHECK(cudaMemsetAsync(out_counts+state.cells,0,(ncell-state.cells)*4,queue));CHECK(cudaMemsetAsync(out_starts+state.cells,0,(ncell-state.cells)*4,queue));}
 if(genesis_pair_prepare_mode()==2){std::vector<double> reference(4ull*atoms);std::vector<int> ids(atoms);
 CHECK(cudaMemcpyAsync(reference.data(),out,reference.size()*8,cudaMemcpyDeviceToHost,queue));CHECK(cudaMemcpyAsync(ids.data(),classes,ids.size()*4,cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));
 if(std::memcmp(reference.data(),coords,reference.size()*8) || std::memcmp(ids.data(),host_cls,ids.size()*4)){std::fprintf(stderr,"GPU pair preparation differs from literal CPU\n");std::abort();}
 std::printf("GENESIS_PAIR_PREPARE: exact translated coordinates charges and classes pass\n");}
 if(!ready)CHECK(cudaEventCreateWithFlags(&ready,cudaEventDisableTiming));CHECK(cudaEventRecord(ready,queue));
 prepared=true;fresh=true;destination=out;class_destination=classes;active=false;return 1;
}

extern "C" int genesis_pair_prepare_force_ready(const double *target,int atoms,int update){
 return update && prepared && fresh && target==destination && atoms==total;
}

extern "C" int genesis_pair_prepare_class_copy(int *target,int atoms,std::uintptr_t stream){
 if(!prepared)return 0;if(atoms!=total)std::abort();auto queue=reinterpret_cast<cudaStream_t>(stream);
 CHECK(cudaStreamWaitEvent(queue,ready,0));CHECK(cudaMemcpyAsync(target,class_destination,4ull*atoms,cudaMemcpyDeviceToDevice,queue));return 1;
}
extern "C" int genesis_pair_prepare_coord_copy(void *target,int atoms,int components,std::uintptr_t stream,int consumer){
 if(!prepared || !fresh || consumer==2)return 0;if(atoms!=total || components<3 || components>4)std::abort();auto queue=reinterpret_cast<cudaStream_t>(stream);
 CHECK(cudaStreamWaitEvent(queue,ready,0));CHECK(cudaMemcpyAsync(target,destination,8ull*components*atoms,cudaMemcpyDeviceToDevice,queue));return 1;
}
extern "C" void genesis_pair_prepare_expire(){fresh=false;}
