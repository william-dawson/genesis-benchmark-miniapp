#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess) { \
  std::fprintf(stderr,"pair sort: %s: %s\n",#call,cudaGetErrorString(e)); std::abort(); } } while(0)

extern "C" void genesis_pair_sort_view(int **,int **,void **);
namespace {
bool flag(const char *name) {
  const char *value=std::getenv(name);
  if(!value || !std::strcmp(value,"0")) return false;
  if(!std::strcmp(value,"1")) return true;
  std::fprintf(stderr,"pair sort: invalid %s\n",name);std::abort();
}
int *keys=nullptr,*sorted_keys=nullptr,*values=nullptr,*order=nullptr,*nonzero=nullptr;
int capacity=0,active=0,checked_nonzero=0;
void *scratch=nullptr;
size_t scratch_capacity=0;
std::vector<int> checked_keys,checked_order;

__global__ void pair_sort_keys(const int *ix,const int *iy,int *key,int *value,
                              int total,int local,int near,int width,int maxload) {
  const int p=blockIdx.x*blockDim.x+threadIdx.x;
  if(p>=total)return;
  const int x=ix[p],y=iy[p];
  int load=16*(((x+7)/8)*((y+3)/4))+(x+7)/8;
  if(p<local)load=maxload-width+(x+y)/2;
  else if(p<near)load=min(load,maxload-width);
  else load=min(load,maxload-width-1);
  key[p]=load;value[p]=p+1;
}
__global__ void pair_sort_count(const int *sorted,int total,int *result) {
  int lo=0,hi=total;
  while(lo<hi) {
    const int mid=lo+(hi-lo)/2;
    if(sorted[mid]>0)lo=mid+1;else hi=mid;
  }
  *result=lo;
}
}

extern "C" int genesis_pair_sort_enabled() {
  static const bool selected=flag("GENESIS_PAIR_SORT");
  return selected;
}
extern "C" int genesis_pair_sort_copy(int *destination,size_t bytes,cudaStream_t stream) {
  if(!genesis_pair_sort_enabled())return 0;
  if(active<=0 || bytes<size_t(active)*sizeof(int))std::abort();
  CHECK(cudaMemcpyAsync(destination,order,size_t(active)*sizeof(int),cudaMemcpyDeviceToDevice,stream));
  return 1;
}
extern "C" void genesis_pair_sort_build(int total,int local,int near,int width,int maxload,
    int *host_nonzero,int *host_ix,int *host_iy,int *host_loads,int *host_order) {
  if(!genesis_pair_sort_enabled() || total<=0 || local<0 || near<local || total<near)std::abort();
  int *ix,*iy;void *queue;
  genesis_pair_sort_view(&ix,&iy,&queue);
  auto stream=reinterpret_cast<cudaStream_t>(queue);
  if(total>capacity) {
    if(keys) {
      CHECK(cudaFree(keys));CHECK(cudaFree(sorted_keys));CHECK(cudaFree(values));CHECK(cudaFree(order));
    }
    capacity=total;
    CHECK(cudaMalloc(&keys,size_t(capacity)*sizeof(int)));
    CHECK(cudaMalloc(&sorted_keys,size_t(capacity)*sizeof(int)));
    CHECK(cudaMalloc(&values,size_t(capacity)*sizeof(int)));
    CHECK(cudaMalloc(&order,size_t(capacity)*sizeof(int)));
  }
  if(!nonzero)CHECK(cudaMalloc(&nonzero,sizeof(int)));
  size_t bytes=0;
  // Stable descending keys retain ascending original IDs within each load bin.
  CHECK(cub::DeviceRadixSort::SortPairsDescending(nullptr,bytes,keys,sorted_keys,values,order,total,0,32,stream));
  if(bytes>scratch_capacity) {
    if(scratch)CHECK(cudaFree(scratch));
    CHECK(cudaMalloc(&scratch,bytes));scratch_capacity=bytes;
  }
  pair_sort_keys<<<(total+255)/256,256,0,stream>>>(ix,iy,keys,values,total,local,near,width,maxload);
  CHECK(cudaGetLastError());
  CHECK(cub::DeviceRadixSort::SortPairsDescending(scratch,bytes,keys,sorted_keys,values,order,total,0,32,stream));
  pair_sort_count<<<1,1,0,stream>>>(sorted_keys,total,nonzero);
  CHECK(cudaGetLastError());
  CHECK(cudaMemcpyAsync(host_nonzero,nonzero,sizeof(int),cudaMemcpyDeviceToHost,stream));
  static const bool verify=flag("GENESIS_PAIR_SORT_VERIFY");
  static const bool poison=flag("GENESIS_PAIR_SORT_POISON");
  if(verify) {
    checked_keys.resize(total);checked_order.resize(total);
    CHECK(cudaMemcpyAsync(host_ix,ix,size_t(total)*sizeof(int),cudaMemcpyDeviceToHost,stream));
    CHECK(cudaMemcpyAsync(host_iy,iy,size_t(total)*sizeof(int),cudaMemcpyDeviceToHost,stream));
    CHECK(cudaMemcpyAsync(checked_keys.data(),keys,size_t(total)*sizeof(int),cudaMemcpyDeviceToHost,stream));
    CHECK(cudaMemcpyAsync(checked_order.data(),order,size_t(total)*sizeof(int),cudaMemcpyDeviceToHost,stream));
  }
  CHECK(cudaStreamSynchronize(stream));
  active=total;checked_nonzero=*host_nonzero;
  if(poison && !verify) {
    std::fill_n(host_ix,total,-123456789);std::fill_n(host_iy,total,-123456789);
    std::fill_n(host_loads,total,-123456789);std::fill_n(host_order,total,-123456789);
    std::puts("GENESIS_PAIR_SORT_POISON: CPU counts loads and order poisoned");
  }
  static bool reported=false;
  if(!reported) {
    std::printf("GENESIS_PAIR_SORT: pairs=%d local=%d near=%d width=%d maxload=%d count_bytes=4\n",
                total,local,near,width,maxload);reported=true;
  }
}
extern "C" void genesis_pair_sort_check(int total,int host_nonzero,const int *loads,const int *indices) {
  if(total!=active || checked_nonzero!=host_nonzero || int(checked_keys.size())!=total ||
     std::memcmp(checked_keys.data(),loads,size_t(total)*sizeof(int)) ||
     std::memcmp(checked_order.data(),indices,size_t(total)*sizeof(int))) {
    std::fputs("GPU pair sorting differs from ordered CPU reference\n",stderr);std::abort();
  }
  std::printf("GENESIS_PAIR_SORT_VERIFY: exact keys count and stable order pass %d %d\n",total,host_nonzero);
}
