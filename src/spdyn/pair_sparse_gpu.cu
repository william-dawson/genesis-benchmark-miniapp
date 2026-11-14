#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess) { \
  fprintf(stderr,"pair sparse: %s: %s\n",#call,cudaGetErrorString(e)); abort(); } } while(0)

static int4 *device_records=nullptr;
static int *device_inverse=nullptr;
static size_t record_capacity=0, inverse_capacity=0;
static int4 *host_records=nullptr;
static size_t host_capacity=0;

extern "C" void *genesis_pair_sparse_staging(int count) {
  const size_t bytes=size_t(count>0 ? count:1)*sizeof(int4);
  if(bytes>host_capacity) {
    if(host_records) CHECK(cudaFreeHost(host_records));
    host_capacity=bytes+bytes/4;
    CHECK(cudaHostAlloc(&host_records,host_capacity,cudaHostAllocPortable));
  }
  return host_records;
}

__device__ __forceinline__ void clear_mask(uint8_t *packed, uint8_t *dense, size_t index) {
  atomicAnd(reinterpret_cast<unsigned *>(packed)+(index>>5),~(1u<<(index&31)));
  atomicAnd(reinterpret_cast<unsigned *>(dense)+(index>>2),~(0xffu<<((index&3)*8)));
}

__global__ void pair_sparse_diagonal(uint8_t *packed, uint8_t *dense,
                                     const int *natom, const int *pairs,
                                     int width, int pitch, int local) {
  const int plane=blockIdx.x;
  if(plane>=local) return;
  const int count=natom[pairs[2*plane]-1];
  for(int atom=threadIdx.x;atom<count;atom+=blockDim.x)
    clear_mask(packed,dense,size_t(plane)*pitch+size_t(atom)*width+atom);
}

__global__ void pair_sparse_clear(uint8_t *packed, uint8_t *dense,
                                  const int4 *records, const int *inverse,
                                  int count, int width, int pitch, int local) {
  const int entry=blockIdx.x*blockDim.x+threadIdx.x;
  if(entry>=count) return;
  const int4 r=records[entry];
  const int plane=r.z<0 ? inverse[-r.z-1]-1 : local+r.z-1;
  const size_t base=size_t(plane)*pitch;
  clear_mask(packed,dense,base+size_t(r.y-1)*width+r.x-1);
  if(r.z<0)
    clear_mask(packed,dense,base+size_t(r.x-1)*width+r.y-1);
}

extern "C" void genesis_pair_sparse_build(uint8_t *packed, uint8_t *dense,
    const int *natom, const int *pairs, const int4 *records, const int *inverse,
    int count, int width, int pitch, int local, int near, cudaStream_t stream) {
  const size_t record_bytes=size_t(count)*sizeof(int4);
  const size_t inverse_bytes=size_t(local)*sizeof(int);
  if(record_bytes>record_capacity) {
    if(device_records) CHECK(cudaFree(device_records));
    CHECK(cudaMalloc(&device_records,record_bytes));
    record_capacity=record_bytes;
  }
  if(inverse_bytes>inverse_capacity) {
    if(device_inverse) CHECK(cudaFree(device_inverse));
    CHECK(cudaMalloc(&device_inverse,inverse_bytes));
    inverse_capacity=inverse_bytes;
  }
  const size_t bytes=size_t(pitch)*near;
  const size_t pack_bytes=(bytes+7)/8;
  CHECK(cudaMemsetAsync(packed,0xff,(pack_bytes+3)&~size_t(3),stream));
  CHECK(cudaMemsetAsync(dense,1,pack_bytes*8,stream));
  if(count) CHECK(cudaMemcpyAsync(device_records,records,record_bytes,cudaMemcpyHostToDevice,stream));
  CHECK(cudaMemcpyAsync(device_inverse,inverse,inverse_bytes,cudaMemcpyHostToDevice,stream));
  if(bytes%8) CHECK(cudaMemsetAsync(packed+pack_bytes-1,(1u<<(bytes%8))-1,1,stream));
  pair_sparse_diagonal<<<local,128,0,stream>>>(packed,dense,natom,pairs,width,pitch,local);
  CHECK(cudaGetLastError());
  if(count) {
    pair_sparse_clear<<<(count+255)/256,256,0,stream>>>(packed,dense,device_records,device_inverse,
                                                   count,width,pitch,local);
    CHECK(cudaGetLastError());
  }
  static int reported=0;
  if(!reported++) fprintf(stdout,"GENESIS_PAIR_SPARSE: records=%d record_bytes=%zu dense_bytes=%zu packed_bytes=%zu\n",
                          count,record_bytes,bytes,pack_bytes);
}
