#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "bonded_incoming_view.h"
extern "C" int genesis_bonded_incoming_cooperative(const MoveView*,int*);
extern "C" int genesis_bonded_hbm_prepare(const MoveView*,int,MoveView*);
extern "C" void genesis_bonded_hbm_publish(const MoveView*,const MoveView*);
namespace {
void require(bool ok,const char *message) {
    if(!ok){std::fprintf(stderr,"GENESIS bonded incoming: %s\n",message);std::abort();}
}
void checked(cudaError_t code){require(code==cudaSuccess,cudaGetErrorString(code));}
int *error_device=nullptr;
struct Snapshot {
    bool pending=false;
    int kind=0;
    MoveView view;
    std::vector<int> counts,candidate_counts;
    std::vector<unsigned char> before[12],candidate[12];
} snapshot;
__device__ void from_packet(MoveField field,int cell,int destination,int source) {
    for(int row=0;row<field.rows;++row) {
        const int out=row+field.rows*(destination+field.capacity*cell);
        const int in=field.source_offset+row+field.source_rows*(source+(field.source_capacity<0?0:field.source_capacity*cell));
        if(field.bytes==8)static_cast<double*>(field.target)[out]=static_cast<const double*>(field.source)[in];
        else if(field.bytes==4)static_cast<int*>(field.target)[out]=static_cast<const int*>(field.source)[in];
        else static_cast<std::int8_t*>(field.target)[out]=static_cast<std::int8_t>(static_cast<const int*>(field.source)[in]);
    }
}
__device__ void from_slot(MoveField field,int cell,int destination,int source) {
    for(int row=0;row<field.rows;++row) {
        const int out=row+field.rows*(destination+field.capacity*cell),in=row+field.rows*(source+field.capacity*cell);
        if(field.bytes==8)static_cast<double*>(field.target)[out]=static_cast<double*>(field.target)[in];
        else if(field.bytes==4)static_cast<int*>(field.target)[out]=static_cast<int*>(field.target)[in];
        else static_cast<std::int8_t*>(field.target)[out]=static_cast<std::int8_t*>(field.target)[in];
    }
}
__global__ void insert_records(MoveView view,int *error) {
    const int cell=blockIdx.x*blockDim.x+threadIdx.x;
    if(cell>=view.cells)return;
    const int old=view.counts[cell],add=view.adds[cell],exits=view.exits[cell];
    if(old<0||add<0||exits<0||exits>old||exits>view.index_capacity){atomicExch(error,1);return;}
    for(int f=0;f<view.fields;++f) {
        const MoveField field=view.field[f];
        if(old>field.capacity||old+add-exits>field.capacity||(field.source_capacity>0 && add>field.source_capacity)){atomicExch(error,1);return;}
    }
    for(int k=0;k<add;++k) {
        const int slot=k<exits?view.indices[k+view.index_capacity*cell]-1:k+old-exits;
        if(slot<0||slot>=(k<exits?old:old+add-exits)){atomicExch(error,1);return;}
        int packet=k;
        if(view.field[11].source)packet+=static_cast<const int*>(view.field[11].source)[2*(cell+view.field[11].source_offset*view.cells)]-1;
        for(int f=0;f<view.fields;++f)from_packet(view.field[f],cell,slot,packet);
    }
    if(add<exits) {
        int removed=0,source=old-1,next=add;
        while(removed<exits-add) {
            bool insert=true;
            for(int k=next;k<exits;++k)if(source==view.indices[k+view.index_capacity*cell]-1){insert=false;++removed;break;}
            if(insert) {
                const int slot=view.indices[next+view.index_capacity*cell]-1;
                if(source<0||slot<0||slot>=old){atomicExch(error,1);return;}
                for(int f=0;f<view.fields;++f)from_slot(view.field[f],cell,slot,source);
                ++removed;++next;
            }
            --source;
        }
    }
    view.counts[cell]=old+add-exits;
}
}
extern "C" int genesis_bonded_incoming_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_BONDED_INCOMING_GPU");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");return int(value);
    }();return mode;
}
extern "C" void genesis_bonded_incoming_start(const MoveView *view,int kind,int verify) {
    require(!snapshot.pending&&view->cells>0&&view->fields>0&&view->fields<=12,"invalid incoming request");
    static const bool supported=[] {
        int pageable=0,ats=0;
        checked(cudaDeviceGetAttribute(&pageable,cudaDevAttrPageableMemoryAccess,0));
        checked(cudaDeviceGetAttribute(&ats,cudaDevAttrPageableMemoryAccessUsesHostPageTables,0));
        require(pageable&&ats,"coherent host page-table access is required");return true;
    }();(void)supported;
    if(verify) {
        snapshot.view=*view;snapshot.kind=kind;
        snapshot.counts.assign(view->counts,view->counts+view->cells);
        for(int f=0;f<view->fields;++f) {
            const auto field=view->field[f];const auto *bytes=static_cast<const unsigned char*>(field.target);
            snapshot.before[f].assign(bytes,bytes+size_t(field.rows)*field.capacity*view->cells*field.bytes);
        }
    }
    if(!error_device)checked(cudaMalloc(&error_device,sizeof(int)));
    checked(cudaMemsetAsync(error_device,0,sizeof(int),0));
    MoveView device=*view;const bool hbm=genesis_bonded_hbm_prepare(view,kind,&device);
    if(!genesis_bonded_incoming_cooperative(&device,error_device))
    insert_records<<<(view->cells+127)/128,128>>>(device,error_device);checked(cudaGetLastError());
    if(hbm)genesis_bonded_hbm_publish(&device,view);
    int error=0;checked(cudaMemcpyAsync(&error,error_device,sizeof(int),cudaMemcpyDeviceToHost,0));checked(cudaStreamSynchronize(0));
    require(!error,"invalid record index/capacity");
    if(verify) {
        snapshot.candidate_counts.assign(view->counts,view->counts+view->cells);
        for(int f=0;f<view->fields;++f) {
            auto field=view->field[f];const auto *bytes=static_cast<const unsigned char*>(field.target);
            snapshot.candidate[f].assign(bytes,bytes+snapshot.before[f].size());
            std::memcpy(field.target,snapshot.before[f].data(),snapshot.before[f].size());
        }
        std::memcpy(view->counts,snapshot.counts.data(),view->cells*sizeof(int));snapshot.pending=true;
    }
}
extern "C" void genesis_bonded_incoming_finish(int kind) {
    if(!snapshot.pending)return;
    require(kind==snapshot.kind,"wrong incoming oracle kind");
    const auto &view=snapshot.view;
    for(int cell=0;cell<view.cells;++cell) {
        const int count=view.counts[cell];require(count==snapshot.candidate_counts[cell],"record count differs from CPU");
        for(int f=0;f<view.fields;++f) {
            const auto field=view.field[f];const size_t stride=size_t(field.rows)*field.capacity*field.bytes;
            const auto *reference=static_cast<const unsigned char*>(field.target)+stride*cell;
            require(!std::memcmp(reference,snapshot.candidate[f].data()+stride*cell,size_t(field.rows)*count*field.bytes),
                "ordered record payload/identity differs from CPU");
        }
    }
    snapshot.pending=false;
    std::printf("GENESIS_BONDED_INCOMING: exact payload and identity pass %d\n",kind);
}
