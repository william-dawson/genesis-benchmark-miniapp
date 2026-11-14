#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "../config.h"
#if defined(_SINGLE) || defined(_MIXED)
#define REAL float
#else
#define REAL double
#endif
#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess) { \
  fprintf(stderr,"pair compact: %s: %s\n",#call,cudaGetErrorString(e)); abort(); } } while(0)

__global__ void pair_compact_build(const REAL *coordinates,const int8_t *moves,
    const int *natom,const int *starts,const int *pairs,uint8_t *xlist,uint8_t *ylist,
    int *nx,int *ny,int atoms,int ld,int cells,int count,REAL radius2,REAL sx,REAL sy,REAL sz) {
  const int plane=blockIdx.x*blockDim.y+threadIdx.y;
  if(plane>=count) return;
  const int lane=threadIdx.x, warp=threadIdx.y;
  const int i=pairs[2*plane]-1, j=pairs[2*plane+1]-1;
  const int ni=natom[i], nj=natom[j], si=starts[i], sj=starts[j];
  __shared__ uint8_t active_x[4][256],active_y[4][256];
  for(int index=lane;index<256;index+=32) {
    active_x[warp][index]=0;active_y[warp][index]=0;
  }
  __syncwarp();
  const size_t move=3*(size_t(j)+size_t(cells)*i);
  const REAL dx=REAL(moves[move])*sx, dy=REAL(moves[move+1])*sy, dz=REAL(moves[move+2])*sz;
  for(int xb=0;xb<ni;xb+=8) {
    const int x=xb+(lane>>2);
    REAL rx=0,ry=0,rz=0;
    if(x<ni) {
      rx=dx+coordinates[si+x];
      ry=dy+coordinates[atoms+si+x];
      rz=dz+coordinates[2*atoms+si+x];
    }
    bool any_x=false;
    for(int yb=0;yb<nj;yb+=4) {
      const int y=yb+(lane&3);
      bool close=false;
      if(x<ni && y<nj) {
        const REAL a=rx-coordinates[sj+y];
        const REAL b=ry-coordinates[atoms+sj+y];
        const REAL c=rz-coordinates[2*atoms+sj+y];
        close=a*a+b*b+c*c<radius2;
      }
      any_x|=close;
      const unsigned hits=__ballot_sync(0xffffffff,close);
      if(lane<4 && yb+lane<nj && (hits&(0x11111111u<<lane))) active_y[warp][yb+lane]=1;
    }
    const unsigned hits=__ballot_sync(0xffffffff,any_x);
    if((lane&3)==0 && x<ni && ((hits>>(lane&~3))&15)) active_x[warp][x]=1;
  }
  __syncwarp();
  int total=0;
  for(int base=0;base<ni;base+=32) {
    const int atom=base+lane;
    const bool keep=atom<ni && active_x[warp][atom];
    const unsigned votes=__ballot_sync(0xffffffff,keep);
    if(keep) xlist[size_t(plane)*ld+total+__popc(votes&((1u<<lane)-1))]=uint8_t(atom+1);
    total+=__popc(votes);
  }
  if(lane==0) nx[plane]=total;
  total=0;
  for(int base=0;base<nj;base+=32) {
    const int atom=base+lane;
    const bool keep=atom<nj && active_y[warp][atom];
    const unsigned votes=__ballot_sync(0xffffffff,keep);
    if(keep) ylist[size_t(plane)*ld+total+__popc(votes&((1u<<lane)-1))]=uint8_t(atom+1);
    total+=__popc(votes);
  }
  if(lane==0) ny[plane]=total;
}

__global__ void pair_compact_check(const uint8_t *a,const uint8_t *b,const int *na,const int *nb,
                                  int ld,int count,int *failures) {
  const int plane=blockIdx.x;
  if(plane>=count) return;
  __shared__ uint8_t aa[256],bb[256];
  const int lane=threadIdx.x;
  for(int k=lane;k<256;k+=blockDim.x) { aa[k]=0;bb[k]=0; }
  __syncthreads();
  if(lane==0 && na[plane]!=nb[plane]) atomicAdd(failures,1);
  for(int k=lane;k<na[plane];k+=blockDim.x) aa[a[size_t(plane)*ld+k]]=1;
  for(int k=lane;k<nb[plane];k+=blockDim.x) bb[b[size_t(plane)*ld+k]]=1;
  __syncthreads();
  for(int k=lane;k<256;k+=blockDim.x) if(aa[k]!=bb[k]) atomicAdd(failures,1);
}

extern "C" void genesis_pair_compact_launch(const REAL *coordinates,const int8_t *moves,
    const int *natom,const int *starts,const int *pairs,uint8_t *xlist,uint8_t *ylist,int *nx,int *ny,
    int atoms,int ld,int local,int boundary,int cells,int count,int capacity,
    REAL radius2,REAL cutoff2,REAL sx,REAL sy,REAL sz,cudaStream_t stream,int verify) {
  static uint8_t *xx=nullptr,*yy=nullptr;
  static int *nnx=nullptr,*nny=nullptr,*failures=nullptr;
  static size_t allocated=0;
  if(verify && size_t(count)*ld>allocated) {
    if(xx) { CHECK(cudaFree(xx));CHECK(cudaFree(yy));CHECK(cudaFree(nnx));CHECK(cudaFree(nny));CHECK(cudaFree(failures)); }
    CHECK(cudaMalloc(&xx,size_t(count)*ld));CHECK(cudaMalloc(&yy,size_t(count)*ld));
    CHECK(cudaMalloc(&nnx,size_t(count)*sizeof(int)));CHECK(cudaMalloc(&nny,size_t(count)*sizeof(int)));
    CHECK(cudaMalloc(&failures,sizeof(int)));allocated=size_t(count)*ld;
  }
  pair_compact_build<<<(count+3)/4,dim3(32,4),0,stream>>>(coordinates,moves,natom,starts,pairs,
      verify?xx:xlist,verify?yy:ylist,verify?nnx:nx,verify?nny:ny,atoms,ld,cells,count,radius2,sx,sy,sz);
  CHECK(cudaGetLastError());
  if(verify) {
    CHECK(cudaMemsetAsync(failures,0,sizeof(int),stream));
    pair_compact_check<<<count,128,0,stream>>>(xlist,xx,nx,nnx,ld,count,failures);CHECK(cudaGetLastError());
    pair_compact_check<<<count,128,0,stream>>>(ylist,yy,ny,nny,ld,count,failures);CHECK(cudaGetLastError());
    int bad=0;
    CHECK(cudaMemcpyAsync(&bad,failures,sizeof(int),cudaMemcpyDeviceToHost,stream));
    CHECK(cudaStreamSynchronize(stream));
    if(bad) { fprintf(stderr,"Pair-list membership verification failed: %d\n",bad);abort(); }
    CHECK(cudaMemcpyAsync(xlist,xx,size_t(count)*ld,cudaMemcpyDeviceToDevice,stream));
    CHECK(cudaMemcpyAsync(ylist,yy,size_t(count)*ld,cudaMemcpyDeviceToDevice,stream));
    CHECK(cudaMemcpyAsync(nx,nnx,size_t(count)*sizeof(int),cudaMemcpyDeviceToDevice,stream));
    CHECK(cudaMemcpyAsync(ny,nny,size_t(count)*sizeof(int),cudaMemcpyDeviceToDevice,stream));
    fprintf(stdout,"GENESIS_PAIR_COMPACT_VERIFY: both atom sets pass\n");
  }
}
