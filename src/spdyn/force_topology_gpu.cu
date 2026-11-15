#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cstdint>
#include "bonded_incoming_view.h"
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"force topology: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_bonded_hbm_working_view(int,MoveView*);
extern "C" int genesis_exclusion_force_view(int,const int**,const int**,const int**,int*);
namespace {
struct Output {int kind,cells,expected,capacity,*offset,*owner,*local,*counts,*starts,*list,*pbc,*period;double *real[4];int8_t *move;};
Output outputs[6]{};int8_t *cell_move=nullptr;int move_pitch=0,*error=nullptr;void *scratch=nullptr;size_t scratch_bytes=0;bool active=false,pending=false;
struct Buffer {void *pointer;size_t bytes;std::vector<unsigned char> candidate;};std::vector<Buffer> buffers;
int mode(){static int m=[](){const char *s=std::getenv("GENESIS_FORCE_TOPOLOGY_GPU");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;if(!std::strcmp(s,"2"))return 2;std::abort();}();return m;}
__global__ void pack_fields(MoveView v,Output out,int *bad){int c=blockIdx.x,n=v.counts[c],base=out.offset[c];
 if(!threadIdx.x && c==v.cells-1 && base+n!=out.expected)atomicExch(bad,1);
 for(int j=threadIdx.x;j<n;j+=blockDim.x){int k=base+j;if(k<0 || k>=out.capacity){atomicExch(bad,2);continue;}
 out.owner[k]=c+1;out.local[k]=j+1;
 auto list=v.field[0];for(int a=0;a<list.rows;++a)out.list[a+list.rows*k]=static_cast<const int*>(list.target)[a+list.rows*(j+list.capacity*c)];
 int pf=out.kind==1?2:(out.kind==3?3:2),rf=out.kind==1?3:(out.kind==3?4:3);auto pbc=v.field[pf];
 for(int a=0;a<pbc.rows;++a)out.pbc[a+pbc.rows*k]=static_cast<const int*>(pbc.target)[a+pbc.rows*(j+pbc.capacity*c)];
 int reals=out.kind==2?4:2;
 for(int a=0;a<reals;++a){auto f=v.field[rf+a];out.real[a][k]=static_cast<const double*>(f.target)[j+f.capacity*c];}
 if(out.kind==3){auto f=v.field[1];out.period[k]=static_cast<const int*>(f.target)[j+f.capacity*c];}
 }
}
__global__ void pack_exclusions(const int *counts,const int *list,int capacity,Output out,const int8_t *move,int pitch,int *bad){int c=blockIdx.x,n=counts[c],base=out.offset[c];
 if(!threadIdx.x && c==out.cells-1 && base+n!=out.expected)atomicExch(bad,3);
 for(int j=threadIdx.x;j<n;j+=blockDim.x){int k=base+j;if(k<0 || k>=out.capacity){atomicExch(bad,4);continue;}
 out.owner[k]=c+1;out.local[k]=j+1;int4 r=reinterpret_cast<const int4*>(list)[j+capacity*c];reinterpret_cast<int4*>(out.list)[k]=r;
 for(int a=0;a<3;++a)out.move[3*k+a]=move[a+3ull*(r.y-1+size_t(pitch)*(r.x-1))];
 }
}
void remember(void *pointer,size_t bytes,cudaStream_t queue){Buffer b{pointer,bytes,{}};b.candidate.resize(bytes);CHECK(cudaMemcpyAsync(b.candidate.data(),pointer,bytes,cudaMemcpyDeviceToHost,queue));buffers.push_back(std::move(b));}
void snapshot(Output o,cudaStream_t q){remember(o.offset,4ull*o.cells,q);remember(o.counts,4ull*o.cells,q);remember(o.owner,4ull*o.expected,q);remember(o.local,4ull*o.expected,q);
 int rows=o.kind==1?2:(o.kind==2?3:4);remember(o.list,4ull*rows*o.expected,q);
 if(o.kind<=4){remember(o.pbc,4ull*(o.kind==1?1:3)*o.expected,q);for(int a=0;a<(o.kind==2?4:2);++a)remember(o.real[a],8ull*o.expected,q);if(o.kind==3)remember(o.period,4ull*o.expected,q);}
 else {remember(o.move,3ull*o.expected,q);remember(o.starts,4ull*o.cells,q);}
}
}
extern "C" void genesis_force_topology_fields(int kind,int cells,int expected,int capacity,int *offset,int *owner,int *local,int *counts,int *list,int *pbc,int *period,double *r0,double *r1,double *r2,double *r3){
 if(kind<1 || kind>4)std::abort();outputs[kind-1]={kind,cells,expected,capacity,offset,owner,local,counts,nullptr,list,pbc,period,{r0,r1,r2,r3},nullptr};
}
extern "C" void genesis_force_topology_exclusion(int kind,int cells,int expected,int capacity,int *offset,int *owner,int *local,int *counts,int *starts,int *list,int8_t *move,const int8_t *host_move,int pitch){
 if(kind<5 || kind>6)std::abort();outputs[kind-1]={kind,cells,expected,capacity,offset,owner,local,counts,starts,list,nullptr,nullptr,{nullptr,nullptr,nullptr,nullptr},move};
 if(!cell_move && mode()){move_pitch=pitch;CHECK(cudaMalloc(&cell_move,3ull*pitch*pitch));CHECK(cudaMemcpy(cell_move,host_move,3ull*pitch*pitch,cudaMemcpyHostToDevice));}if(cell_move && pitch!=move_pitch)std::abort();
}
extern "C" void genesis_force_topology_start(std::uintptr_t stream){
 active=false;pending=false;buffers.clear();if(!mode())return;MoveView views[4];for(int k=0;k<4;++k)if(!genesis_bonded_hbm_working_view(k+1,&views[k]))return;
 const int *nc[2],*lists[2],*starts[2];int caps[2];for(int k=0;k<2;++k)if(!genesis_exclusion_force_view(k,&nc[k],&lists[k],&starts[k],&caps[k]))return;
 auto queue=reinterpret_cast<cudaStream_t>(stream);if(!error)CHECK(cudaMalloc(&error,4));CHECK(cudaMemsetAsync(error,0,4,queue));
 for(int k=0;k<6;++k){auto o=outputs[k];if(o.kind!=k+1 || o.cells<1 || o.expected<0 || o.expected>o.capacity)std::abort();
 const int *count=k<4?views[k].counts:nc[k-4];size_t bytes=0;CHECK(cub::DeviceScan::ExclusiveSum(nullptr,bytes,count,o.offset,o.cells,queue));
 if(bytes>scratch_bytes){if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_bytes=bytes;}
 CHECK(cub::DeviceScan::ExclusiveSum(scratch,bytes,count,o.offset,o.cells,queue));CHECK(cudaMemcpyAsync(o.counts,count,4ull*o.cells,cudaMemcpyDeviceToDevice,queue));
 if(k<4)pack_fields<<<o.cells,128,0,queue>>>(views[k],o,error);
 else {CHECK(cudaMemcpyAsync(o.starts,starts[k-4],4ull*o.cells,cudaMemcpyDeviceToDevice,queue));pack_exclusions<<<o.cells,128,0,queue>>>(count,lists[k-4],caps[k-4],o,cell_move,move_pitch,error);}
 CHECK(cudaGetLastError());if(mode()==2)snapshot(o,queue);
 }
 int bad=0;CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));if(bad){std::fprintf(stderr,"GPU force topology count/capacity error %d\n",bad);std::abort();}
 active=true;pending=mode()==2;
}
extern "C" int genesis_force_topology_bypass(){return active && mode()==1;}
extern "C" void genesis_force_topology_check(std::uintptr_t stream){if(!pending)return;auto queue=reinterpret_cast<cudaStream_t>(stream);
 for(const auto &b:buffers){std::vector<unsigned char> reference(b.bytes);CHECK(cudaMemcpyAsync(reference.data(),b.pointer,b.bytes,cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));if(std::memcmp(reference.data(),b.candidate.data(),b.bytes)){std::fprintf(stderr,"GPU force topology differs from original CPU packing\n");std::abort();}}
 pending=false;std::printf("GENESIS_FORCE_TOPOLOGY: exact compact indices counts and parameters pass\n");
}
