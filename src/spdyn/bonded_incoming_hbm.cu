#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "bonded_incoming_view.h"
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"bonded incoming HBM: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_improper_packets(MoveView*);
extern "C" void genesis_bonded_device_packets(const int**,const int**,const int**,const int**,const int**,const double**);
namespace {
struct Cache {bool ready=false;MoveView v{},host{};};Cache cache[5];
__global__ void packet_metadata(MoveView v,const int *src,const int *dst,const int *slots,int kind,int *error){int c=blockIdx.x,b=c+(kind-1)*v.cells;
 int start=src[2*b]-1,n=start<0?0:src[2*b+1]-start,first=dst[2*b]-1,m=first<0?0:dst[2*b+1]-first;
 if(n>v.index_capacity || n<0 || m<0){if(!threadIdx.x)atomicExch(error,1);return;}
 if(!threadIdx.x){const_cast<int*>(v.exits)[c]=n;const_cast<int*>(v.adds)[c]=m;}
 for(int j=threadIdx.x;j<n;j+=blockDim.x)const_cast<int*>(v.indices)[j+v.index_capacity*c]=slots[start+j];
}
__global__ void publish(MoveView d,MoveView h){int c=blockIdx.x,n=d.counts[c];if(!threadIdx.x)h.counts[c]=n;
 for(int f=0;f<d.fields;++f){auto a=d.field[f],b=h.field[f];if(n<0 || n>a.capacity)return;
 for(int i=threadIdx.x;i<n*a.rows*a.bytes;i+=blockDim.x)static_cast<unsigned char*>(b.target)[size_t(c)*b.rows*b.capacity*b.bytes+i]=static_cast<const unsigned char*>(a.target)[size_t(c)*a.rows*a.capacity*a.bytes+i];}
}
int mode(){static int result=[](){const char *s=std::getenv("GENESIS_BONDED_INCOMING_HBM");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;std::abort();}();return result;}
void *clone(const void *p,size_t bytes){void *out;CHECK(cudaMalloc(&out,bytes?bytes:1));if(bytes)CHECK(cudaMemcpyAsync(out,p,bytes,cudaMemcpyHostToDevice,0));return out;}
}
extern "C" int genesis_bonded_hbm_prepare(const MoveView *host,int kind,MoveView *out){
 if(!mode())return 0;if(kind<1 || kind>5 || host->fields>=12)std::abort();auto &c=cache[kind-1];c.host=*host;
 if(!c.ready){c.v=*host;c.v.counts=static_cast<int*>(clone(host->counts,4ull*host->cells));
 for(int f=0;f<host->fields;++f){auto a=host->field[f];c.v.field[f].target=clone(a.target,size_t(a.rows)*a.capacity*a.bytes*host->cells);}
 c.v.adds=static_cast<const int*>(clone(host->adds,4ull*host->cells));c.v.exits=static_cast<const int*>(clone(host->exits,4ull*host->cells));c.v.indices=static_cast<const int*>(clone(host->indices,4ull*host->index_capacity*host->cells));c.ready=true;}
 if(c.v.cells!=host->cells || c.v.fields!=host->fields)std::abort();MoveView v=c.v;
 if(kind<=3){const char *resident=std::getenv("GENESIS_BONDED_RESIDENCY");if(!resident || std::strcmp(resident,"1"))std::abort();
 const int *src,*dst,*slots,*ids,*iv;const double *rv;genesis_bonded_device_packets(&src,&dst,&slots,&ids,&iv,&rv);
 static int *error=nullptr;if(!error)CHECK(cudaMalloc(&error,4));CHECK(cudaMemsetAsync(error,0,4,0));
 packet_metadata<<<v.cells,128>>>(v,src,dst,slots,kind,error);CHECK(cudaGetLastError());int bad=0;CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,0));CHECK(cudaStreamSynchronize(0));if(bad)std::abort();
 for(int f=0;f<v.fields;++f){auto &a=v.field[f];a.source=f==v.fields-1?static_cast<const void*>(ids):(a.bytes==8?static_cast<const void*>(rv):static_cast<const void*>(iv));a.source_rows=f==v.fields-1?1:(a.bytes==8?4:9);a.source_capacity=-1;}
 v.field[11].source=dst;v.field[11].source_offset=kind-1;
 }else if(kind==4 && genesis_improper_packets(&v)){
 }else{
 CHECK(cudaMemcpyAsync(const_cast<int*>(v.adds),host->adds,4ull*v.cells,cudaMemcpyHostToDevice,0));CHECK(cudaMemcpyAsync(const_cast<int*>(v.exits),host->exits,4ull*v.cells,cudaMemcpyHostToDevice,0));CHECK(cudaMemcpyAsync(const_cast<int*>(v.indices),host->indices,4ull*v.cells*v.index_capacity,cudaMemcpyHostToDevice,0));
 for(int f=0;f<v.fields;++f){auto h=host->field[f];int shared=-1;
 for(int j=0;j<f;++j)if(host->field[j].source==h.source){shared=j;break;}
 if(shared>=0)c.v.field[f].source=c.v.field[shared].source;
 else {size_t bytes=size_t(h.source_rows)*h.source_capacity*v.cells*(h.bytes==8?8:4);
 if(!c.v.field[f].source || c.v.field[f].source==h.source)c.v.field[f].source=clone(h.source,bytes);
 else CHECK(cudaMemcpyAsync(const_cast<void*>(c.v.field[f].source),h.source,bytes,cudaMemcpyHostToDevice,0));}
 v.field[f].source=c.v.field[f].source;}
 }
 *out=v;return 1;
}
extern "C" void genesis_bonded_hbm_publish(const MoveView *device,const MoveView *host){
 static int mirror=[](){const char *s=std::getenv("GENESIS_BONDED_HOST_MIRROR");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;if(!std::strcmp(s,"2"))return 2;std::abort();}();
 bool retire=mirror!=0;
 for(const char *name:{"GENESIS_BONDED_INCOMING_GPU","GENESIS_BONDED_INCOMING_HBM","GENESIS_BONDED_RESIDENCY","GENESIS_EXCLUSION_GPU","GENESIS_FORCE_TOPOLOGY_GPU","GENESIS_IMPROPER_OUTGOING_GPU"}){const char *v=std::getenv(name);if(!v || std::strcmp(v,"1"))retire=false;}
 for(const char *name:{"GENESIS_BONDED_MIGRATION_VERIFY","GENESIS_EXCLUSION_VERIFY","GENESIS_FORCE_VERIFY","GENESIS_BONDED_VERIFY"}){const char *v=std::getenv(name);if(v && std::strcmp(v,"0"))retire=false;}
 for(int k=0;k<4;++k)if(!cache[k].ready)retire=false;
 if(!retire){publish<<<host->cells,128>>>(*device,*host);CHECK(cudaGetLastError());return;}
 if(mirror==2){for(int f=0;f<host->fields;++f){auto a=host->field[f];std::memset(a.target,0xff,size_t(a.rows)*a.capacity*a.bytes*host->cells);}std::memset(host->counts,0xff,4ull*host->cells);std::printf("GENESIS_BONDED_HOST_MIRROR: CPU fields and counts poisoned\n");}
}
extern "C" int genesis_bonded_hbm_metadata(int kind,const int **counts,const int **ids){if(!mode() || kind<1 || kind>3 || !cache[kind-1].ready)return 0;auto &c=cache[kind-1];*counts=c.v.counts;*ids=static_cast<const int*>(c.v.field[c.v.fields-1].target);return 1;}
extern "C" int genesis_bonded_hbm_exclusion(int kind,const int **counts,const int **list,const signed char **types,const int **period){
 if(!mode() || kind<1 || kind>3 || !cache[kind-1].ready)return 0;auto &v=cache[kind-1].v;
 *counts=v.counts;*list=static_cast<const int*>(v.field[0].target);*types=static_cast<const signed char*>(v.field[kind==3?2:1].target);*period=kind==3?static_cast<const int*>(v.field[1].target):nullptr;return 1;
}
extern "C" int genesis_bonded_hbm_working_view(int kind,MoveView *out){if(!mode() || kind<1 || kind>5 || !cache[kind-1].ready)return 0;*out=cache[kind-1].v;return 1;}

extern "C" int genesis_bonded_hbm_host_view(int kind,MoveView *out){if(!mode() || kind<1 || kind>5 || !cache[kind-1].ready)return 0;*out=cache[kind-1].host;return 1;}

extern "C" int genesis_bonded_hbm_ready(int kind){return mode() && kind>=1 && kind<=5 && cache[kind-1].ready;}
extern "C" void genesis_bonded_hbm_initialize(const MoveView *host,int kind){MoveView device;if(!genesis_bonded_hbm_prepare(host,kind,&device))std::abort();}
