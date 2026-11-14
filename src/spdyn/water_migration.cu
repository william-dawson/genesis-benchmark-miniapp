#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include "../config.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <chrono>

extern "C" void genesis_water_scatter(double *,double *,int *,int *,int *,int *,
    const double *,const int *,const int *,int,int,int);

namespace {
void checked(cudaError_t code) {
    if(code!=cudaSuccess) {std::fprintf(stderr,"Water migration CUDA: %s\n",cudaGetErrorString(code));std::abort();}
}
void require(bool condition,const char *message) {
    if(!condition) {std::fprintf(stderr,"Water migration: %s\n",message);std::abort();}
}
int flag(const char *name) {
    const char *v=std::getenv(name);
    if(!v)return 0;
    require(!std::strcmp(v,"0")||!std::strcmp(v,"1"),name);
    return v[0]=='1';
}
struct State {
    const double *x=nullptr,*v=nullptr;
    const int *water=nullptr,*counts=nullptr;
    int width=0,cells=0,water_width=0,capacity=0,globals=0,total=0;
    int *keys[2]{},*values[2]{},*prefix=nullptr,*ids=nullptr,*ranges=nullptr,*error=nullptr,*packed_ids=nullptr;
    std::int16_t *map=nullptr;
    double *packed=nullptr,*host=nullptr;
    int *host_ids=nullptr,*host_ranges=nullptr;
    void *temporary=nullptr;
    size_t temporary_bytes=0;
    cudaEvent_t ready=nullptr,start=nullptr,finish=nullptr;
    bool latched=false,timing=false;
} state;

__device__ int cell_coordinate(double value,double origin,double box,double cell,int count) {
    double shifted=__dadd_rn(value,-origin);
    const double offset=__dadd_rn(__dmul_rn(box,0.5),-__dmul_rn(box,round(shifted/box)));
    shifted=__dadd_rn(shifted,offset);
    int index=int(shifted/cell);
    if(index==count)--index;
    return index;
}
__global__ void classify(const double *x,const int *water,const int *counts,const int *prefix,
    const std::int16_t *map,int *keys,int *values,int *error,int width,int cells,int water_width,
    double3 box,double3 origin,double3 size,int3 grid) {
    int entry=blockIdx.x*blockDim.x+threadIdx.x;
    const int cell=entry/water_width,water_index=entry%water_width;
    if(cell>=cells||water_index>=counts[cell])return;
    int atom=water[3ull*entry]-1;
    if(atom<0||atom>=width) {atomicExch(error,1);return;}
    const double *position=x+3ull*(atom+width*cell);
    const int cx=cell_coordinate(position[0],origin.x,box.x,size.x,grid.x);
    const int cy=cell_coordinate(position[1],origin.y,box.y,size.y,grid.y);
    const int cz=cell_coordinate(position[2],origin.z,box.z,size.z,grid.z);
    const int output=prefix[cell]+water_index;
    values[output]=entry;keys[output]=0;
    if(cx<0||cx>=grid.x||cy<0||cy>=grid.y||cz<0||cz>=grid.z) {atomicExch(error,2);return;}
    const int destination=int(map[cx+grid.x*(cy+grid.y*cz)])-1;
    if(destination<0||destination>=cells) {atomicExch(error,3);return;}
    keys[output]=destination==cell?cell:cells+destination;
}

__global__ void pack(const double *x,const double *v,const int *water,const int *ids,
    const int *keys,const int *values,int *ranges,double *packed,int *packed_ids,int *error,
    int width,int water_width,int total) {
    const int item=blockIdx.x*blockDim.x+threadIdx.x;
    if(item>=total)return;
    const int key=keys[item];
    if(item==0||keys[item-1]!=key)ranges[2*key]=item;
    if(item==total-1||keys[item+1]!=key)ranges[2*key+1]=item+1;
    const int entry=values[item],cell=entry/water_width;
    for(int list=0;list<3;++list) {
        const int atom=water[3ull*entry+list]-1;
        if(atom<0||atom>=width) {atomicExch(error,4);continue;}
        const size_t input=atom+size_t(width)*cell;
        packed_ids[3ull*item+list]=ids[input];
        for(int component=0;component<3;++component) {
            packed[18ull*item+6*list+component]=x[3*input+component];
            packed[18ull*item+6*list+3+component]=v[3*input+component];
        }
    }
}
}

extern "C" int genesis_water_enabled() {
    static int enabled=flag("GENESIS_WATER_GPU");
#if defined(_SINGLE) || defined(_MIXED)
    require(!enabled,"requires double precision");
#endif
    return enabled;
}
extern "C" int genesis_water_verify() {static int enabled=flag("GENESIS_WATER_VERIFY");return enabled;}
extern "C" void genesis_water_latch(const double *x,const double *v,const int *water,const int *counts,
    int width,int cells,int water_width) {
    if(!genesis_water_enabled())return;
    auto &s=state;
    if(!s.ready)checked(cudaEventCreateWithFlags(&s.ready,cudaEventDisableTiming));
    if(s.width)require(s.width==width&&s.cells==cells&&s.water_width==water_width,"shape changed");
    s.x=x;s.v=v;s.water=water;s.counts=counts;s.width=width;s.cells=cells;s.water_width=water_width;
    checked(cudaEventRecord(s.ready,0));s.latched=true;
}

extern "C" int genesis_water_outgoing(const int *counts,const int *ids,const std::int16_t *map,
    const double *box,const double *origin,const double *size,const int *grid,int *stay,int *move,
    double *stay_real,double *move_real,int *stay_ids,int *move_ids,int stay_capacity,int move_capacity) {
    auto &s=state;
    if(!genesis_water_enabled()||!s.latched)return 0;
    s.latched=false;
    const auto begin=std::chrono::steady_clock::now();
    const int globals=grid[0]*grid[1]*grid[2];
    std::vector<int> prefix(s.cells+1);
    for(int c=0;c<s.cells;++c) {
        require(counts[c]>=0&&counts[c]<=s.water_width,"invalid water count");
        prefix[c+1]=prefix[c]+counts[c];
    }
    s.total=prefix.back();
    require(s.total>0&&globals>0,"empty shape");
    if(!s.capacity) {
        s.capacity=s.water_width*s.cells;s.globals=globals;
        for(auto &p:s.keys)checked(cudaMalloc(&p,size_t(s.capacity)*sizeof(int)));
        for(auto &p:s.values)checked(cudaMalloc(&p,size_t(s.capacity)*sizeof(int)));
        checked(cudaMalloc(&s.prefix,size_t(s.cells+1)*sizeof(int)));
        checked(cudaMalloc(&s.ids,size_t(s.width)*s.cells*sizeof(int)));
        checked(cudaMalloc(&s.map,size_t(globals)*sizeof(std::int16_t)));
        checked(cudaMalloc(&s.ranges,4ull*s.cells*sizeof(int)));
        checked(cudaMalloc(&s.error,sizeof(int)));
        checked(cudaMalloc(&s.packed,18ull*s.capacity*sizeof(double)));
        checked(cudaMalloc(&s.packed_ids,3ull*s.capacity*sizeof(int)));
        checked(cudaMallocHost(&s.host,18ull*s.capacity*sizeof(double)));
        checked(cudaMallocHost(&s.host_ids,3ull*s.capacity*sizeof(int)));
        checked(cudaMallocHost(&s.host_ranges,4ull*s.cells*sizeof(int)));
        checked(cub::DeviceRadixSort::SortPairs(nullptr,s.temporary_bytes,s.keys[0],s.keys[1],
            s.values[0],s.values[1],s.capacity));
        checked(cudaMalloc(&s.temporary,s.temporary_bytes));
        s.timing=flag("GENESIS_WATER_TIMING");
        if(s.timing) {checked(cudaEventCreate(&s.start));checked(cudaEventCreate(&s.finish));}
    }
    require(globals==s.globals,"global cell shape changed");
    checked(cudaStreamWaitEvent(0,s.ready,0));
    checked(cudaMemcpy(s.prefix,prefix.data(),size_t(s.cells+1)*sizeof(int),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.ids,ids,size_t(s.width)*s.cells*sizeof(int),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.map,map,size_t(globals)*sizeof(std::int16_t),cudaMemcpyHostToDevice));
    checked(cudaMemset(s.ranges,0xff,4ull*s.cells*sizeof(int)));
    checked(cudaMemset(s.error,0,sizeof(int)));
    if(s.timing)checked(cudaEventRecord(s.start,0));
    classify<<<(s.capacity+255)/256,256>>>(s.x,s.water,s.counts,s.prefix,s.map,s.keys[0],s.values[0],s.error,
        s.width,s.cells,s.water_width,make_double3(box[0],box[1],box[2]),
        make_double3(origin[0],origin[1],origin[2]),make_double3(size[0],size[1],size[2]),
        make_int3(grid[0],grid[1],grid[2]));
    checked(cudaGetLastError());
    int bits=0;while((1ull<<bits)<2ull*s.cells)++bits;
    checked(cub::DeviceRadixSort::SortPairs(s.temporary,s.temporary_bytes,s.keys[0],s.keys[1],
        s.values[0],s.values[1],s.total,0,bits));
    pack<<<(s.total+255)/256,256>>>(s.x,s.v,s.water,s.ids,s.keys[1],s.values[1],s.ranges,s.packed,
        s.packed_ids,s.error,s.width,s.water_width,s.total);
    checked(cudaGetLastError());
    if(s.timing)checked(cudaEventRecord(s.finish,0));
    checked(cudaMemcpy(s.host,s.packed,18ull*s.total*sizeof(double),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(s.host_ids,s.packed_ids,3ull*s.total*sizeof(int),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(s.host_ranges,s.ranges,4ull*s.cells*sizeof(int),cudaMemcpyDeviceToHost));
    int error=0;checked(cudaMemcpy(&error,s.error,sizeof(int),cudaMemcpyDeviceToHost));
    require(!error,"invalid classification or atom index");
    const auto copied=std::chrono::steady_clock::now();
    for(int bucket=0;bucket<2*s.cells;++bucket) {
        const int first=s.host_ranges[2*bucket],last=s.host_ranges[2*bucket+1];
        const int count=first<0?0:last-first;
        const int capacity=bucket<s.cells?stay_capacity:move_capacity;
        require(count>=0&&count<=capacity,"destination buffer overflow");
    }
    genesis_water_scatter(stay_real,move_real,stay_ids,move_ids,stay,move,s.host,s.host_ids,s.host_ranges,
        s.cells,stay_capacity,move_capacity);
    if(s.timing) {
        float ms=0;checked(cudaEventElapsedTime(&ms,s.start,s.finish));
        const auto end=std::chrono::steady_clock::now();
        std::printf("GENESIS_WATER_TIMING: %d %.9f %.9f %.9f\n",s.total,double(ms)*1e-3,
            std::chrono::duration<double>(copied-begin).count(),std::chrono::duration<double>(end-copied).count());
    }
    return 1;
}

extern "C" void genesis_water_check(const int *stay,const int *move,const double *stay_real,const double *move_real,
    const int *stay_ids,const int *move_ids,int stay_capacity,int move_capacity) {
    auto &s=state;
    int seen=0;
    for(int bucket=0;bucket<2*s.cells;++bucket) {
        const int first=s.host_ranges[2*bucket],last=s.host_ranges[2*bucket+1];
        const int count=first<0?0:last-first,cell=bucket%s.cells;
        const int capacity=bucket<s.cells?stay_capacity:move_capacity;
        require((bucket<s.cells?stay:move)[cell]==count,"CPU/GPU water counts differ");
        if(!count)continue;
        const double *real=(bucket<s.cells?stay_real:move_real)+18ull*capacity*cell;
        const int *integer=(bucket<s.cells?stay_ids:move_ids)+3ull*capacity*cell;
        require(!std::memcmp(real,s.host+18ull*first,18ull*count*sizeof(double)),"CPU/GPU water coordinates or velocities differ");
        require(!std::memcmp(integer,s.host_ids+3ull*first,3ull*count*sizeof(int)),"CPU/GPU water identity/order differs");
        seen+=count;
    }
    require(seen==s.total,"water conservation failed");
    std::printf("GENESIS_WATER_VERIFY: exact packing pass %d\n",seen);
}
