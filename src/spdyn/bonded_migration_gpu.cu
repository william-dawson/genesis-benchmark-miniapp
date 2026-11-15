#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess) { \
  std::fprintf(stderr,"bonded migration: %s: %s\n",#call,cudaGetErrorString(e));std::abort(); } } while(0)
namespace {
struct Cache { int *integers;double *reals;int total,sites; };
struct CacheSet { Cache kind[3]; } caches{};
bool ready[3]={false,false,false},initialized=false;
int4 *input=nullptr,*selected=nullptr,*host_input=nullptr;
uint8_t *flags=nullptr;
int *result=nullptr,*host_result=nullptr,*keys=nullptr,*sorted_keys=nullptr,*values=nullptr,*order=nullptr;
int *exits=nullptr,*ids=nullptr,*integers=nullptr,*source_ranges=nullptr,*destination_ranges=nullptr;
double *reals=nullptr;
int *host_exits=nullptr,*host_ids=nullptr,*host_integers=nullptr,*host_source_ranges=nullptr,*host_destination_ranges=nullptr;
double *host_reals=nullptr;
int16_t *mapping=nullptr,*owners=nullptr;
int capacity=0,cached_atoms=0,cached_cells=0,cached_pitch=0;
const int16_t *owner_address=nullptr;
void *scratch=nullptr,*sort_scratch=nullptr;
size_t scratch_capacity=0,sort_capacity=0;
cudaStream_t queue=nullptr;
cudaEvent_t begin_event,copy_event,end_event;
bool flag(const char *name) {
  const char *v=std::getenv(name);
  if(!v || !std::strcmp(v,"0"))return false;
  if(!std::strcmp(v,"1"))return true;
  std::fprintf(stderr,"Invalid %s\n",name);std::abort();
}
__global__ void classify_bonded(int4 *records,uint8_t *active,const int16_t *map,const int16_t *owner,
                              int n,int atoms,int cells,int pitch,int *error,CacheSet data) {
  const int j=blockIdx.x*blockDim.x+threadIdx.x;
  if(j>=n)return;
  const int4 r=records[j];const int kind=r.z>>16,source=r.z&65535,id=r.x-1;
  if(kind<0 || kind>2 || source<1 || source>cells) {
    atomicExch(error,1);active[j]=0;return;
  }
  Cache c;
  c.integers=kind==0?data.kind[0].integers:(kind==1?data.kind[1].integers:data.kind[2].integers);
  c.total=kind==0?data.kind[0].total:(kind==1?data.kind[1].total:data.kind[2].total);
  c.sites=kind==0?data.kind[0].sites:(kind==1?data.kind[1].sites:data.kind[2].sites);
  if(id<0 || id>=c.total) { atomicExch(error,1);active[j]=0;return; }
  const int a=c.integers[9*id],b=c.integers[9*id+c.sites-1];
  if(a<1 || a>atoms || b<1 || b>atoms) { atomicExch(error,2);active[j]=0;return; }
  const int first=map[2*(a-1)],last=map[2*(b-1)];
  if(first<1 || first>cells || last<1 || last>cells) { atomicExch(error,3);active[j]=0;return; }
  const int destination=owner[first-1+size_t(pitch)*(last-1)];
  if(destination<1 || destination>cells) { atomicExch(error,4);active[j]=0;return; }
  records[j]=make_int4(r.z,r.w,destination,r.x);
  active[j]=destination!=source;
}
__global__ void bonded_sort_keys(const int4 *moves,int *key,int *value,int n,int cells) {
  const int j=blockIdx.x*blockDim.x+threadIdx.x;
  if(j>=n)return;
  key[j]=(moves[j].x>>16)*cells+moves[j].z-1;value[j]=j;
}
__global__ void gather_bonded(const int4 *moves,const int *sorted,const int *key,int n,int cells,
    int *out_exits,int *out_ids,int *out_int,double *out_real,int *src,int *dst,CacheSet data) {
  const int j=blockIdx.x*blockDim.x+threadIdx.x;
  if(j>=n)return;
  const int4 source=moves[j];const int s=(source.x>>16)*cells+(source.x&65535)-1;
  out_exits[j]=source.y;
  if(j==0 || moves[j-1].x!=source.x)src[2*s]=j+1;
  if(j==n-1 || moves[j+1].x!=source.x)src[2*s+1]=j+1;
  const int4 r=moves[sorted[j]];const int kind=r.x>>16,id=r.w-1;
  Cache c;
  c.integers=kind==0?data.kind[0].integers:(kind==1?data.kind[1].integers:data.kind[2].integers);
  c.reals=kind==0?data.kind[0].reals:(kind==1?data.kind[1].reals:data.kind[2].reals);
  out_ids[j]=r.w;
  #pragma unroll
  for(int f=0;f<9;++f)out_int[9*j+f]=c.integers[9*id+f];
  #pragma unroll
  for(int f=0;f<4;++f)out_real[4*j+f]=c.reals[4*id+f];
  const int k=key[j];
  if(j==0 || key[j-1]!=k)dst[2*k]=j+1;
  if(j==n-1 || key[j+1]!=k)dst[2*k+1]=j+1;
}
}
extern "C" void *genesis_bonded_staging(int n) {
  if(!initialized) {
    CHECK(cudaStreamCreateWithFlags(&queue,cudaStreamNonBlocking));
    CHECK(cudaMalloc(&result,2*sizeof(int)));CHECK(cudaHostAlloc(&host_result,2*sizeof(int),cudaHostAllocPortable));
    CHECK(cudaEventCreate(&begin_event));CHECK(cudaEventCreate(&copy_event));CHECK(cudaEventCreate(&end_event));
    initialized=true;
  }
  if(n>capacity) {
    if(input)std::abort();
    capacity=n+n/4+1;
    CHECK(cudaMalloc(&input,size_t(capacity)*sizeof(int4)));CHECK(cudaMalloc(&selected,size_t(capacity)*sizeof(int4)));
    CHECK(cudaMalloc(&flags,size_t(capacity)));CHECK(cudaHostAlloc(&host_input,size_t(capacity)*sizeof(int4),cudaHostAllocPortable));
    CHECK(cudaMalloc(&keys,size_t(capacity)*4));CHECK(cudaMalloc(&sorted_keys,size_t(capacity)*4));
    CHECK(cudaMalloc(&values,size_t(capacity)*4));CHECK(cudaMalloc(&order,size_t(capacity)*4));
    CHECK(cudaMalloc(&exits,size_t(capacity)*4));CHECK(cudaMalloc(&ids,size_t(capacity)*4));
    CHECK(cudaMalloc(&integers,size_t(capacity)*9*4));CHECK(cudaMalloc(&reals,size_t(capacity)*4*sizeof(double)));
    CHECK(cudaHostAlloc(&host_exits,size_t(capacity)*4,cudaHostAllocPortable));
    CHECK(cudaHostAlloc(&host_ids,size_t(capacity)*4,cudaHostAllocPortable));
    CHECK(cudaHostAlloc(&host_integers,size_t(capacity)*9*4,cudaHostAllocPortable));
    CHECK(cudaHostAlloc(&host_reals,size_t(capacity)*4*sizeof(double),cudaHostAllocPortable));
  }
  return host_input;
}
extern "C" void genesis_bonded_cache(int kind,int n,const int *integer_data,const double *real_data) {
  if(!initialized || kind<1 || kind>3 || ready[kind-1] || n<0)std::abort();
  Cache &c=caches.kind[kind-1];c.total=n;c.sites=kind+1;
  CHECK(cudaMalloc(&c.integers,size_t(n?n:1)*9*4));CHECK(cudaMalloc(&c.reals,size_t(n?n:1)*4*sizeof(double)));
  if(n) {
    CHECK(cudaMemcpyAsync(c.integers,integer_data,size_t(n)*9*4,cudaMemcpyHostToDevice,queue));
    CHECK(cudaMemcpyAsync(c.reals,real_data,size_t(n)*4*sizeof(double),cudaMemcpyHostToDevice,queue));
  }
  ready[kind-1]=true;
}
extern "C" void genesis_bonded_batch(int n,const int16_t *map,const int16_t *owner,int atoms,int cells,int pitch,
    int *moved,void **out_exits,void **out_ids,void **out_integers,void **out_reals,void **out_src,void **out_dst) {
  const auto start=std::chrono::steady_clock::now();static const bool timing=flag("GENESIS_BONDED_MIGRATION_TIMING");
  if(!initialized || n<0 || n>capacity || cells<=0 || cells>=32768 || pitch<cells)std::abort();
  int total=0;for(int k=0;k<3;++k) { if(!ready[k])std::abort();total+=caches.kind[k].total; }
  if(n!=total)std::abort();
  if(!mapping) {
    CHECK(cudaMalloc(&mapping,size_t(atoms)*2*sizeof(int16_t)));CHECK(cudaMalloc(&owners,size_t(pitch)*pitch*sizeof(int16_t)));
    CHECK(cudaMalloc(&source_ranges,size_t(cells)*6*4));CHECK(cudaMalloc(&destination_ranges,size_t(cells)*6*4));
    CHECK(cudaHostAlloc(&host_source_ranges,size_t(cells)*6*4,cudaHostAllocPortable));
    CHECK(cudaHostAlloc(&host_destination_ranges,size_t(cells)*6*4,cudaHostAllocPortable));
    cached_atoms=atoms;cached_cells=cells;cached_pitch=pitch;owner_address=owner;
    CHECK(cudaMemcpyAsync(owners,owner,size_t(pitch)*pitch*sizeof(int16_t),cudaMemcpyHostToDevice,queue));
  }
  if(atoms!=cached_atoms || cells!=cached_cells || pitch!=cached_pitch || owner!=owner_address)std::abort();
  if(timing)CHECK(cudaEventRecord(begin_event,queue));
  CHECK(cudaMemcpyAsync(mapping,map,size_t(atoms)*2*sizeof(int16_t),cudaMemcpyHostToDevice,queue));
  if(n)CHECK(cudaMemcpyAsync(input,host_input,size_t(n)*sizeof(int4),cudaMemcpyHostToDevice,queue));
  CHECK(cudaMemsetAsync(result,0,2*sizeof(int),queue));
  if(timing)CHECK(cudaEventRecord(copy_event,queue));
  if(n) {
    size_t bytes=0;CHECK(cub::DeviceSelect::Flagged(nullptr,bytes,input,flags,selected,result,n,queue));
    if(bytes>scratch_capacity) { if(scratch)CHECK(cudaFree(scratch));CHECK(cudaMalloc(&scratch,bytes));scratch_capacity=bytes; }
    classify_bonded<<<(n+255)/256,256,0,queue>>>(input,flags,mapping,owners,n,atoms,cells,pitch,result+1,caches);
    CHECK(cudaGetLastError());CHECK(cub::DeviceSelect::Flagged(scratch,bytes,input,flags,selected,result,n,queue));
  }
  CHECK(cudaMemcpyAsync(host_result,result,2*sizeof(int),cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));
  if(host_result[1] || host_result[0]<0 || host_result[0]>n) {
    std::fprintf(stderr,"Invalid bonded batch: error=%d moved=%d\n",host_result[1],host_result[0]);std::abort();
  }
  *moved=host_result[0];const int m=*moved;
  CHECK(cudaMemsetAsync(source_ranges,0,size_t(cells)*6*4,queue));
  CHECK(cudaMemsetAsync(destination_ranges,0,size_t(cells)*6*4,queue));
  if(m) {
    bonded_sort_keys<<<(m+255)/256,256,0,queue>>>(selected,keys,values,m,cells);CHECK(cudaGetLastError());
    int bits=0;for(unsigned k=3*cells-1;k;k>>=1)++bits;
    size_t bytes=0;CHECK(cub::DeviceRadixSort::SortPairs(nullptr,bytes,keys,sorted_keys,values,order,m,0,bits,queue));
    if(bytes>sort_capacity) { if(sort_scratch)CHECK(cudaFree(sort_scratch));CHECK(cudaMalloc(&sort_scratch,bytes));sort_capacity=bytes; }
    CHECK(cub::DeviceRadixSort::SortPairs(sort_scratch,bytes,keys,sorted_keys,values,order,m,0,bits,queue));
    gather_bonded<<<(m+255)/256,256,0,queue>>>(selected,order,sorted_keys,m,cells,exits,ids,integers,reals,
                                          source_ranges,destination_ranges,caches);CHECK(cudaGetLastError());
  }
  if(timing)CHECK(cudaEventRecord(end_event,queue));
  if(m) {
    CHECK(cudaMemcpyAsync(host_exits,exits,size_t(m)*4,cudaMemcpyDeviceToHost,queue));
    CHECK(cudaMemcpyAsync(host_ids,ids,size_t(m)*4,cudaMemcpyDeviceToHost,queue));
    CHECK(cudaMemcpyAsync(host_integers,integers,size_t(m)*9*4,cudaMemcpyDeviceToHost,queue));
    CHECK(cudaMemcpyAsync(host_reals,reals,size_t(m)*4*sizeof(double),cudaMemcpyDeviceToHost,queue));
  }
  CHECK(cudaMemcpyAsync(host_source_ranges,source_ranges,size_t(cells)*6*4,cudaMemcpyDeviceToHost,queue));
  CHECK(cudaMemcpyAsync(host_destination_ranges,destination_ranges,size_t(cells)*6*4,cudaMemcpyDeviceToHost,queue));
  CHECK(cudaStreamSynchronize(queue));
  *out_exits=host_exits;*out_ids=host_ids;*out_integers=host_integers;*out_reals=host_reals;
  *out_src=host_source_ranges;*out_dst=host_destination_ranges;
  if(timing) {
    float upload,pipeline;CHECK(cudaEventElapsedTime(&upload,begin_event,copy_event));
    CHECK(cudaEventElapsedTime(&pipeline,copy_event,end_event));
    const double host=std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-start).count();
    std::printf("GENESIS_BONDED_MIGRATION_GPU: 0 %d %d %.6f %.6f %.6f\n",n,m,upload*1000,pipeline*1000,host);
  }
}
