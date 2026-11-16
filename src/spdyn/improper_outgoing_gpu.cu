#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "bonded_incoming_view.h"
#include "migration_views.h"
using genesis::WorkingReal;
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"improper outgoing: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_bonded_hbm_working_view(int,MoveView*);
extern "C" int genesis_bonded_hbm_host_view(int,MoveView*);
extern "C" int genesis_incoming_pair_view(ParticleView*,const genesis::WorkingReal**);
namespace {
MoveView device{},host{};bool packets=false,pending=false;int *prefix=nullptr,*result=nullptr,*error=nullptr;int16_t *owners=nullptr;
int4 *records=nullptr,*selected=nullptr;uint8_t *flags=nullptr;void *scratch=nullptr;size_t scratch_bytes=0;int total=0;
std::vector<int> saved_add,saved_exit,saved_index,saved_int;std::vector<WorkingReal> saved_real;
__global__ void finish_prefix(const int *counts,const int *prefix,int cells,int total,int *bad){if(!threadIdx.x && !blockIdx.x && counts[cells-1]+prefix[cells-1]!=total)atomicExch(bad,1);}
__global__ void classify(MoveView v,const int16_t *map,const int16_t *owner,int pitch,const int *prefix,int4 *out,uint8_t *flags,int *bad){int c=blockIdx.x,n=v.counts[c];
 const int *list=static_cast<const int*>(v.field[0].target);
 for(int j=threadIdx.x;j<n;j+=blockDim.x){int a=list[4*(j+v.field[0].capacity*c)]-1,b=list[4*(j+v.field[0].capacity*c)+3]-1;
 int first=map[2*a]-1,last=map[2*b]-1,dest=-1;
 if(first<0 || first>=v.cells || last<0 || last>=v.cells)atomicExch(bad,2);else dest=owner[first+size_t(pitch)*last]-1;
 if(dest<0 || dest>=v.cells)atomicExch(bad,3);out[prefix[c]+j]=make_int4(c,j,dest,0);flags[prefix[c]+j]=dest!=c;}
 __syncthreads();if(threadIdx.x)return;int nout=0;for(int j=0;j<n;++j)if(flags[prefix[c]+j]){if(nout>=v.index_capacity){atomicExch(bad,4);break;}const_cast<int*>(v.indices)[nout+++v.index_capacity*c]=j+1;}const_cast<int*>(v.exits)[c]=nout;
}
__global__ void pack(MoveView v,const int4 *moves,int moved,int *bad){int c=blockIdx.x;if(threadIdx.x)return;int n=0;
 for(int j=0;j<moved;++j){int4 r=moves[j];if(r.z!=c)continue;
 for(int f=0;f<v.fields;++f){MoveField a=v.field[f];if(n>=a.source_capacity){atomicExch(bad,5);continue;}
 for(int row=0;row<a.rows;++row){int input=row+a.rows*(r.y+a.capacity*r.x),output=a.source_offset+row+a.source_rows*(n+a.source_capacity*c);
 if(a.bytes==8)const_cast<double*>(static_cast<const double*>(a.source))[output]=static_cast<const double*>(a.target)[input];
 else const_cast<int*>(static_cast<const int*>(a.source))[output]=static_cast<const int*>(a.target)[input];}}
 ++n;}const_cast<int*>(v.adds)[c]=n;
}
void publish(){CHECK(cudaMemcpyAsync(const_cast<int*>(host.adds),device.adds,4ull*device.cells,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(const_cast<int*>(host.exits),device.exits,4ull*device.cells,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(const_cast<int*>(host.indices),device.indices,4ull*device.cells*device.index_capacity,cudaMemcpyDeviceToHost,0));
 for(int f:{0,3}){auto a=device.field[f],h=host.field[f];size_t bytes=size_t(a.source_rows)*a.source_capacity*device.cells*(a.bytes==8?8:4);CHECK(cudaMemcpyAsync(const_cast<void*>(h.source),a.source,bytes,cudaMemcpyDeviceToHost,0));}CHECK(cudaStreamSynchronize(0));
 saved_add.assign(host.adds,host.adds+host.cells);saved_exit.assign(host.exits,host.exits+host.cells);saved_index.assign(host.indices,host.indices+size_t(host.cells)*host.index_capacity);
 auto i=host.field[0],r=host.field[3];saved_int.assign(static_cast<const int*>(i.source),static_cast<const int*>(i.source)+size_t(i.source_rows)*i.source_capacity*host.cells);saved_real.assign(static_cast<const WorkingReal*>(r.source),static_cast<const WorkingReal*>(r.source)+size_t(r.source_rows)*r.source_capacity*host.cells);pending=true;
}
}
extern "C" int genesis_improper_outgoing_mode(){static int mode=[](){const char *s=std::getenv("GENESIS_IMPROPER_OUTGOING_GPU");if(!s || !std::strcmp(s,"0"))return 0;if(std::strcmp(s,"1") && std::strcmp(s,"2"))std::abort();const char *v=std::getenv("GENESIS_BONDED_MIGRATION_VERIFY");return s[0]=='2' || (v && !std::strcmp(v,"1"))?2:1;}();return mode;}
extern "C" int genesis_improper_packets(MoveView *view){if(!packets)return 0;*view=device;packets=false;return 1;}
extern "C" int genesis_improper_outgoing(const int16_t *owner,int pitch,int n,int verify){
 if(!genesis_bonded_hbm_working_view(4,&device))return 0;if(!genesis_bonded_hbm_host_view(4,&host))std::abort();ParticleView particle;const genesis::WorkingReal *translation;
 if(!genesis_incoming_pair_view(&particle,&translation))std::abort();if(n<=0 || particle.cells!=device.cells || pitch<device.cells)std::abort();
 if(!prefix){total=n;CHECK(cudaMalloc(&prefix,4ull*device.cells));CHECK(cudaMalloc(&result,4));CHECK(cudaMalloc(&error,4));CHECK(cudaMalloc(&owners,2ull*pitch*pitch));CHECK(cudaMemcpyAsync(owners,owner,2ull*pitch*pitch,cudaMemcpyHostToDevice,0));CHECK(cudaMalloc(&records,size_t(total)*sizeof(int4)));CHECK(cudaMalloc(&selected,size_t(total)*sizeof(int4)));CHECK(cudaMalloc(&flags,total));}
 if(n!=total)std::abort();CHECK(cudaMemsetAsync(error,0,4,0));size_t scan=0,select=0;CHECK(cub::DeviceScan::ExclusiveSum(nullptr,scan,device.counts,prefix,device.cells));CHECK(cub::DeviceSelect::Flagged(nullptr,select,records,flags,selected,result,total));size_t bytes=scan>select?scan:select;
 if(bytes>scratch_bytes){if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_bytes=bytes;}
 CHECK(cub::DeviceScan::ExclusiveSum(scratch,scan,device.counts,prefix,device.cells));finish_prefix<<<1,1>>>(device.counts,prefix,device.cells,total,error);CHECK(cudaGetLastError());
 classify<<<device.cells,128>>>(device,particle.reverse,owners,pitch,prefix,records,flags,error);CHECK(cudaGetLastError());CHECK(cub::DeviceSelect::Flagged(scratch,select,records,flags,selected,result,total));int moved=0,bad=0;
 CHECK(cudaMemcpyAsync(&moved,result,4,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,0));CHECK(cudaStreamSynchronize(0));if(bad || moved<0 || moved>total)std::abort();
 pack<<<device.cells,32>>>(device,selected,moved,error);CHECK(cudaGetLastError());CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,0));CHECK(cudaStreamSynchronize(0));if(bad)std::abort();packets=true;if(verify)publish();return 1;
}
extern "C" void genesis_improper_outgoing_check(){if(!pending)std::abort();auto iv=host.field[0],rv=host.field[3];
 for(int c=0;c<host.cells;++c){int a=host.adds[c],e=host.exits[c];size_t ib=size_t(iv.source_rows)*iv.source_capacity*c,rb=size_t(rv.source_rows)*rv.source_capacity*c;
 if(a!=saved_add[c] || e!=saved_exit[c] || std::memcmp(host.indices+host.index_capacity*c,saved_index.data()+host.index_capacity*c,size_t(e)*4) || std::memcmp(static_cast<const int*>(iv.source)+ib,saved_int.data()+ib,size_t(a)*iv.source_rows*4) || std::memcmp(static_cast<const WorkingReal*>(rv.source)+rb,saved_real.data()+rb,size_t(a)*rv.source_rows*sizeof(WorkingReal))){std::fprintf(stderr,"Improper outgoing differs from CPU cell %d\n",c);std::abort();}}
 pending=false;std::printf("GENESIS_IMPROPER_OUTGOING: exact counts exits and payload pass\n");
}
