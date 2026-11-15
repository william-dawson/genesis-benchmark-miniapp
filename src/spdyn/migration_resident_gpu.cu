#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "migration_views.h"
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"migration residency: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_incoming_pair_view(ParticleView*,const double**);
namespace {void *scratch=nullptr;size_t scratch_bytes=0;
__global__ void last_prefix(const int *counts,int *prefix,int n){if(!threadIdx.x && !blockIdx.x)prefix[n]=prefix[n-1]+counts[n-1];}
}
extern "C" int genesis_migration_resident(){static int mode=[](){const char *s=std::getenv("GENESIS_MIGRATION_RESIDENT");if(!s || !std::strcmp(s,"0"))return 0;if(std::strcmp(s,"1"))std::abort();const char *gpu=std::getenv("GENESIS_INCOMING_GPU"),*hbm=std::getenv("GENESIS_INCOMING_STORAGE");if(!gpu || gpu[0]=='0' || !hbm || std::strcmp(hbm,"1"))std::abort();return 1;}();return mode;}
extern "C" int genesis_migration_host_packets(){if(!genesis_migration_resident())return 1;
 for(const char *name:{"GENESIS_WATER_VERIFY","GENESIS_HGR_VERIFY","GENESIS_MIGRATION_VERIFY"}){const char *s=std::getenv(name);if(s && !std::strcmp(s,"1"))return 1;}
 const char *mode=std::getenv("GENESIS_INCOMING_GPU");return mode && !std::strcmp(mode,"2");}
extern "C" int genesis_migration_prefix(int kind,int *prefix,int *ids,int cells,int width,int connect){
 if(!genesis_migration_resident())return 0;ParticleView view;const double *translation=nullptr;
 if(!genesis_incoming_pair_view(&view,&translation))return 0;
 if(view.cells!=cells || view.width!=width || (kind==2 && view.connect!=connect))std::abort();
 const int *counts=kind==1?view.nw:view.hlocal;int n=kind==1?cells:cells*connect;
 size_t bytes=0;CHECK(cub::DeviceScan::ExclusiveSum(nullptr,bytes,counts,prefix,n));
 if(bytes>scratch_bytes){if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_bytes=bytes;}
 CHECK(cub::DeviceScan::ExclusiveSum(scratch,bytes,counts,prefix,n));last_prefix<<<1,1>>>(counts,prefix,n);CHECK(cudaGetLastError());
 CHECK(cudaMemcpyAsync(ids,view.ids,size_t(width)*cells*4,cudaMemcpyDeviceToDevice,0));return 1;
}
