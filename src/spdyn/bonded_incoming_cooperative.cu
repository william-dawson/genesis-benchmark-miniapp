#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "bonded_incoming_view.h"
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"cooperative bonded incoming: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
namespace {
__device__ void packet_row(MoveField f,int cell,int to,int from,int row){int dst=row+f.rows*(to+f.capacity*cell),src=f.source_offset+row+f.source_rows*(from+(f.source_capacity<0?0:f.source_capacity*cell));
 if(f.bytes==8)static_cast<double*>(f.target)[dst]=static_cast<const double*>(f.source)[src];
 else if(f.bytes==4)static_cast<int*>(f.target)[dst]=static_cast<const int*>(f.source)[src];
 else static_cast<signed char*>(f.target)[dst]=static_cast<signed char>(static_cast<const int*>(f.source)[src]);
}
__device__ void slot_row(MoveField f,int cell,int to,int from,int row){int dst=row+f.rows*(to+f.capacity*cell),src=row+f.rows*(from+f.capacity*cell);
 if(f.bytes==8)static_cast<double*>(f.target)[dst]=static_cast<const double*>(f.target)[src];
 else if(f.bytes==4)static_cast<int*>(f.target)[dst]=static_cast<const int*>(f.target)[src];
 else static_cast<signed char*>(f.target)[dst]=static_cast<const signed char*>(f.target)[src];
}
__global__ void insert_cell(MoveView v,int *error){int cell=blockIdx.x,old=v.counts[cell],add=v.adds[cell],exits=v.exits[cell],next_count=old+add-exits;
 if(old<0 || add<0 || exits<0 || exits>old || exits>v.index_capacity){if(!threadIdx.x)atomicExch(error,1);return;}
 for(int f=0;f<v.fields;++f){auto a=v.field[f];if(old>a.capacity || next_count>a.capacity || (a.source_capacity>0 && add>a.source_capacity)){if(!threadIdx.x)atomicExch(error,1);return;}}
 int packet_base=0;if(v.field[11].source)packet_base=static_cast<const int*>(v.field[11].source)[2*(cell+v.field[11].source_offset*v.cells)]-1;
 for(int f=0;f<v.fields;++f){auto a=v.field[f];for(int index=threadIdx.x;index<add*a.rows;index+=blockDim.x){int k=index/a.rows,row=index%a.rows,dst=k<exits?v.indices[k+v.index_capacity*cell]-1:k+old-exits;
 if(dst<0 || dst>=(k<exits?old:next_count)){atomicExch(error,1);continue;}packet_row(a,cell,dst,packet_base+k,row);}}
 __syncthreads();extern __shared__ int moves[];__shared__ int count,failed;
 if(!threadIdx.x){count=0;failed=0;
 if(add<exits){int removed=0,source=old-1,next=add;
 while(removed<exits-add){bool insert=true;for(int k=next;k<exits;++k)if(source==v.indices[k+v.index_capacity*cell]-1){insert=false;++removed;break;}
 if(insert){int dst=v.indices[next+v.index_capacity*cell]-1;
 // Tail sources and retained destinations are disjoint after packet insertion.
 if(source<next_count || source>=old || dst<0 || dst>=next_count || count>=v.index_capacity){failed=1;atomicExch(error,2);break;}
 moves[2*count]=source;moves[2*count+1]=dst;++count;++removed;++next;}--source;}
 }}
 __syncthreads();if(failed)return;
 for(int f=0;f<v.fields;++f){auto a=v.field[f];for(int index=threadIdx.x;index<count*a.rows;index+=blockDim.x){int k=index/a.rows,row=index%a.rows;slot_row(a,cell,moves[2*k+1],moves[2*k],row);}}
 __syncthreads();if(!threadIdx.x)v.counts[cell]=next_count;
}
}
extern "C" int genesis_bonded_incoming_cooperative(const MoveView *v,int *error){static int mode=[](){const char *s=std::getenv("GENESIS_BONDED_INCOMING_COOPERATIVE");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;std::abort();}();if(!mode)return 0;
 insert_cell<<<v->cells,128,8ull*v->index_capacity>>>(*v,error);CHECK(cudaGetLastError());return 1;
}
