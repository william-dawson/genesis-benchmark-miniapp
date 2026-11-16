#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "migration_views.h"
#include "migration_coordinates.h"
using genesis::IntegrationReal;
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"atom outgoing: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_incoming_cached_view(const ParticleView*,ParticleView*);
namespace {
AtomRecords device{},host{};ParticleView atoms{};int *prefix=nullptr,*error=nullptr,*result=nullptr;
int4 *records=nullptr,*selected=nullptr;uint8_t *flags=nullptr;int16_t *mapping=nullptr;
int capacity=0,cells=0;void *scratch=nullptr;size_t scratch_bytes=0;bool pending=false;
std::vector<int> saved_add,saved_exit,saved_indices,saved_int;std::vector<IntegrationReal> saved_real;
struct Geometry {double3 box,origin,size;int3 grid;};
template<class T>T* writable(const T *p){return const_cast<T*>(p);}
__device__ int coordinate(double x,double origin,double box,double size,int n){
#if defined(_SINGLE)
 return genesis::single_cell_coordinate(x,origin,box,size,n);
#else
 double shift=__dadd_rn(x,-origin);shift=__dadd_rn(shift,__dadd_rn(__dmul_rn(box,0.5),-__dmul_rn(box,round(shift/box))));int c=int(shift/size);if(c==n)--c;return c;
#endif
}
__global__ void starts(ParticleView v,int *prefix,int *bad,int total){if(threadIdx.x || blockIdx.x)return;int sum=0;
 for(int c=0;c<v.cells;++c){prefix[c]=sum;int n=v.noh[c];if(n<0 || n>v.width)atomicExch(bad,1);sum+=n;}if(sum!=total)atomicExch(bad,2);}
__global__ void classify(ParticleView v,AtomRecords out,const int *prefix,const int16_t *map,Geometry g,int4 *records,uint8_t *flags,int *bad){
 int c=blockIdx.x,n=v.noh[c];
 for(int a=threadIdx.x;a<n;a+=blockDim.x){int index=a+v.width*c;IntegrationReal *x=v.x+3*index;
 int cx=coordinate(x[0],g.origin.x,g.box.x,g.size.x,g.grid.x),cy=coordinate(x[1],g.origin.y,g.box.y,g.size.y,g.grid.y),cz=coordinate(x[2],g.origin.z,g.box.z,g.size.z,g.grid.z);
 int dest=-1;if(cx<0 || cx>=g.grid.x || cy<0 || cy>=g.grid.y || cz<0 || cz>=g.grid.z)atomicExch(bad,3);else dest=map[cx+g.grid.x*(cy+g.grid.y*cz)]-1;
 if(dest<0 || dest>=v.cells)atomicExch(bad,4);
 records[prefix[c]+a]=make_int4(c,a,dest,0);flags[prefix[c]+a]=dest!=c;
 }
 __syncthreads();if(threadIdx.x)return;int exits=0;
 for(int a=0;a<n;++a)if(flags[prefix[c]+a]){if(exits>=out.index_capacity){atomicExch(bad,5);break;}const_cast<int*>(out.indices)[exits+++out.index_capacity*c]=a+1;}
 const_cast<int*>(out.exits)[c]=exits;
}
__global__ void pack(ParticleView v,AtomRecords out,const int4 *moves,int moved,int *bad){int c=blockIdx.x;if(threadIdx.x)return;int add=0;
 for(int j=0;j<moved;++j){int4 r=moves[j];if(r.z!=c)continue;if(add>=out.real_capacity || add>=out.int_capacity){atomicExch(bad,6);break;}
 int source=r.y+v.width*r.x,base=out.real_rows*(add+out.real_capacity*c),ib=out.int_rows*(add+out.int_capacity*c);
 for(int a=0;a<3;++a){const_cast<IntegrationReal*>(out.reals)[base+a]=v.x[3*source+a];const_cast<IntegrationReal*>(out.reals)[base+3+a]=v.v[3*source+a];}
 const_cast<IntegrationReal*>(out.reals)[base+6]=v.q[source];const_cast<IntegrationReal*>(out.reals)[base+7]=v.m[source];
 const_cast<int*>(out.integers)[ib]=v.classes[source];const_cast<int*>(out.integers)[ib+1]=v.ids[source];const_cast<int*>(out.integers)[ib+2]=v.solute[source];++add;
 }const_cast<int*>(out.adds)[c]=add;
}
void snapshot(){saved_add.assign(host.adds,host.adds+cells);saved_exit.assign(host.exits,host.exits+cells);
 saved_indices.assign(host.indices,host.indices+size_t(cells)*host.index_capacity);
 saved_int.assign(host.integers,host.integers+size_t(cells)*host.int_capacity*host.int_rows);
 saved_real.assign(host.reals,host.reals+size_t(cells)*host.real_capacity*host.real_rows);pending=true;}
}
extern "C" int genesis_atom_outgoing_mode(){static int mode=[](){const char *v=std::getenv("GENESIS_ATOM_OUTGOING_GPU");if(!v || !std::strcmp(v,"0"))return 0;if(!std::strcmp(v,"1"))return 1;if(!std::strcmp(v,"2"))return 2;std::abort();}();return mode;}
extern "C" int genesis_atom_packets(AtomRecords *out){if(!genesis_atom_outgoing_mode() || !device.reals)return 0;*out=device;return 1;}
extern "C" void genesis_atom_outgoing(const ParticleView *view,const AtomRecords *output,const double *box,const double *origin,const double *size,const int *grid,const int16_t *map,int verify){
 if(!genesis_incoming_cached_view(view,&atoms)){std::fprintf(stderr,"Atom outgoing requires HBM incoming storage\n");std::abort();}
 host=*output;
 if(!device.reals){cells=view->cells;for(int c=0;c<cells;++c)capacity+=view->noh[c];if(capacity<=0)capacity=0;
 device=host;IntegrationReal *rp;int *ip,*ap,*ep,*xp;
 CHECK(cudaMalloc(&rp,size_t(cells)*host.real_capacity*host.real_rows*sizeof(IntegrationReal)));device.reals=rp;
 CHECK(cudaMalloc(&ip,size_t(cells)*host.int_capacity*host.int_rows*4));device.integers=ip;
 CHECK(cudaMalloc(&ap,cells*4));device.adds=ap;CHECK(cudaMalloc(&ep,cells*4));device.exits=ep;
 CHECK(cudaMalloc(&xp,size_t(cells)*host.index_capacity*4));device.indices=xp;
 CHECK(cudaMalloc(&prefix,cells*4));CHECK(cudaMalloc(&error,4));CHECK(cudaMalloc(&result,4));
 CHECK(cudaMalloc(&records,size_t(capacity?capacity:1)*sizeof(int4)));CHECK(cudaMalloc(&selected,size_t(capacity?capacity:1)*sizeof(int4)));CHECK(cudaMalloc(&flags,capacity?capacity:1));
 size_t globals=size_t(grid[0])*grid[1]*grid[2];CHECK(cudaMalloc(&mapping,globals*2));CHECK(cudaMemcpyAsync(mapping,map,globals*2,cudaMemcpyHostToDevice,0));
 }
 if(cells!=view->cells)std::abort();CHECK(cudaMemsetAsync(error,0,4,0));
 starts<<<1,1>>>(atoms,prefix,error,capacity);CHECK(cudaGetLastError());Geometry g{make_double3(box[0],box[1],box[2]),make_double3(origin[0],origin[1],origin[2]),make_double3(size[0],size[1],size[2]),make_int3(grid[0],grid[1],grid[2])};
 classify<<<cells,128>>>(atoms,device,prefix,mapping,g,records,flags,error);CHECK(cudaGetLastError());
 size_t bytes=0;CHECK(cub::DeviceSelect::Flagged(nullptr,bytes,records,flags,selected,result,capacity));
 if(bytes>scratch_bytes){if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_bytes=bytes;}
 CHECK(cub::DeviceSelect::Flagged(scratch,bytes,records,flags,selected,result,capacity));
 int moved=0,bad=0;CHECK(cudaMemcpyAsync(&moved,result,4,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,0));CHECK(cudaStreamSynchronize(0));
 if(bad || moved<0 || moved>capacity){std::fprintf(stderr,"Atom outgoing classification error %d\n",bad);std::abort();}
 pack<<<cells,32>>>(atoms,device,selected,moved,error);CHECK(cudaGetLastError());
 CHECK(cudaMemcpyAsync(writable(host.adds),device.adds,cells*4,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(writable(host.exits),device.exits,cells*4,cudaMemcpyDeviceToHost,0));
 if(verify){CHECK(cudaMemcpyAsync(writable(host.indices),device.indices,size_t(cells)*host.index_capacity*4,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(writable(host.integers),device.integers,size_t(cells)*host.int_capacity*host.int_rows*4,cudaMemcpyDeviceToHost,0));CHECK(cudaMemcpyAsync(writable(host.reals),device.reals,size_t(cells)*host.real_capacity*host.real_rows*sizeof(IntegrationReal),cudaMemcpyDeviceToHost,0));}
 CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,0));CHECK(cudaStreamSynchronize(0));if(bad)std::abort();if(verify)snapshot();
}
extern "C" void genesis_atom_outgoing_check(){if(!pending)std::abort();
 for(int c=0;c<cells;++c){int a=host.adds[c],e=host.exits[c];if(a!=saved_add[c] || e!=saved_exit[c] || std::memcmp(host.indices+host.index_capacity*c,saved_indices.data()+host.index_capacity*c,size_t(e)*4) || std::memcmp(host.integers+size_t(host.int_rows)*host.int_capacity*c,saved_int.data()+size_t(host.int_rows)*host.int_capacity*c,size_t(a)*host.int_rows*4) || std::memcmp(host.reals+size_t(host.real_rows)*host.real_capacity*c,saved_real.data()+size_t(host.real_rows)*host.real_capacity*c,size_t(a)*host.real_rows*sizeof(IntegrationReal))){std::fprintf(stderr,"Atom outgoing differs from literal CPU at cell %d\n",c);std::abort();}}
 pending=false;std::printf("GENESIS_ATOM_OUTGOING: exact counts exits and payload pass\n");
}
