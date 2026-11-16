#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "precision_types.h"
using genesis::WorkingReal;
#define CHECK(c) do { cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"bonded residency: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_bonded_hbm_metadata(int,const int**,const int**);
namespace {
struct View {int *ids,*counts,*adds,*exits,*indices,*integers;WorkingReal *reals;int cells,idcap,exitcap,ni,intcap,nr,realcap;};
struct Views {View v[3];};
Views views{};
int *counts=nullptr,*offsets=nullptr,*error=nullptr,capacity=0;
void *scratch=nullptr;size_t scratch_bytes=0;
__global__ void record_counts(Views all,int *count,int cells) {
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<3*cells)count[i]=all.v[i/cells].counts[i%cells];
}
__global__ void records(Views all,const int *offset,int4 *out,int cells,int total,int *bad) {
  int b=blockIdx.x,k=b/cells,c=b%cells;View v=all.v[k];int n=v.counts[c];
  if(n<0 || n>v.idcap || offset[b]+n>total || (b==3*cells-1 && offset[b]+n!=total)) {if(threadIdx.x==0)atomicExch(bad,1);return;}
  for(int s=threadIdx.x;s<n;s+=blockDim.x)out[offset[b]+s]=make_int4(v.ids[s+v.idcap*c],0,c+1+(k<<16),s+1);
}
__global__ void scatter(View v,int kind,const int *src,const int *dst,const int *exit_slots,
                        const int *ids,const int *integers,const WorkingReal *reals,int *incoming,int incoming_cap,int *bad) {
  int c=blockIdx.x,b=c+(kind-1)*v.cells;
  int first=src[2*b]-1,n=first>=0?src[2*b+1]-first:0;
  int start=dst[2*b]-1,m=start>=0?dst[2*b+1]-start:0;
  if(n>v.exitcap || m>v.intcap || m>v.realcap || m>incoming_cap){if(threadIdx.x==0)atomicExch(bad,2);return;}
  if(threadIdx.x==0){v.exits[c]=n;v.adds[c]=m;}
  for(int j=threadIdx.x;j<n;j+=blockDim.x)v.indices[j+v.exitcap*c]=exit_slots[first+j];
  for(int j=threadIdx.x;j<m;j+=blockDim.x){
    for(int f=0;f<v.ni;++f)v.integers[f+v.ni*(j+v.intcap*c)]=integers[f+9*(start+j)];
    for(int f=0;f<v.nr;++f)v.reals[f+v.nr*(j+v.realcap*c)]=reals[f+4*(start+j)];
    incoming[j+incoming_cap*c]=ids[start+j];
  }
}
void check(cudaStream_t queue){int bad=0;CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));if(bad){std::fprintf(stderr,"bonded residency bounds %d\n",bad);std::abort();}}
}
extern "C" int genesis_bonded_residency_enabled(){
  static int mode=[](){const char *v=std::getenv("GENESIS_BONDED_RESIDENCY");if(!v || !std::strcmp(v,"0"))return 0;if(!std::strcmp(v,"1"))return 1;std::abort();}();return mode;
}
extern "C" void genesis_bonded_residency_register(int kind,int *ids,int *count,int *adds,int *exits,int *indices,
 int *integers,WorkingReal *reals,int cells,int idcap,int exitcap,int ni,int intcap,int nr,int realcap){
 if(kind<1 || kind>3 || cells<=0)std::abort();
 views.v[kind-1]={ids,count,adds,exits,indices,integers,reals,cells,idcap,exitcap,ni,intcap,nr,realcap};
}
extern "C" int genesis_bonded_residency_prepare(int4 *output,int total,int cells,cudaStream_t queue){
 if(!genesis_bonded_residency_enabled())return 0;
 if(!error)CHECK(cudaMalloc(&error,4));CHECK(cudaMemsetAsync(error,0,4,queue));
 if(!counts){capacity=3*cells;CHECK(cudaMalloc(&counts,4*capacity));CHECK(cudaMalloc(&offsets,4*capacity));}
 if(capacity!=3*cells)std::abort();
 Views current=views;for(int k=0;k<3;++k){const int *count,*ids;if(genesis_bonded_hbm_metadata(k+1,&count,&ids)){current.v[k].counts=const_cast<int*>(count);current.v[k].ids=const_cast<int*>(ids);}}
 record_counts<<<(capacity+255)/256,256,0,queue>>>(current,counts,cells);CHECK(cudaGetLastError());
 size_t bytes=0;CHECK(cub::DeviceScan::ExclusiveSum(nullptr,bytes,counts,offsets,capacity,queue));
 if(bytes>scratch_bytes){if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_bytes=bytes;}
 CHECK(cub::DeviceScan::ExclusiveSum(scratch,bytes,counts,offsets,capacity,queue));
 records<<<capacity,128,0,queue>>>(current,offsets,output,cells,total,error);CHECK(cudaGetLastError());check(queue);return 1;
}
extern "C" void genesis_bonded_residency_scatter(int kind,const int *src,const int *dst,const int *exits,
 const int *ids,const int *integers,const WorkingReal *reals,int *incoming,int cap,cudaStream_t queue){
 CHECK(cudaMemsetAsync(error,0,4,queue));View v=views.v[kind-1];
 scatter<<<v.cells,128,0,queue>>>(v,kind,src,dst,exits,ids,integers,reals,incoming,cap,error);
 CHECK(cudaGetLastError());check(queue);
}
