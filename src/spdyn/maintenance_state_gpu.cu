#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "precision_types.h"
using genesis::IntegrationReal;

namespace {
struct State {
    const IntegrationReal *reference=nullptr;
    int width=0,cells=0;
    bool pending=false;
    std::vector<IntegrationReal> before,after;
} state;
void require(bool condition,const char *message) {
    if(!condition) {std::fprintf(stderr,"GENESIS maintenance state: %s\n",message);std::abort();}
}
void checked(cudaError_t code) {
    if(code!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS maintenance state: %s\n",cudaGetErrorString(code));std::abort();
    }
}
}
extern "C" int genesis_maint_state_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_MAINT_STATE");
        if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");
        return int(value);
    }();
    return mode;
}
extern "C" int genesis_maint_state_begin(const IntegrationReal *reference,int width,int cells,int eligible) {
    if(!genesis_maint_state_mode()||!eligible)return 0;
    auto &s=state;
    require(!s.pending&&reference&&width>0&&cells>0,"invalid or nested maintenance handoff");
    s.reference=reference;s.width=width;s.cells=cells;s.pending=true;
    if(genesis_maint_state_mode()==2) {
        s.before.resize(3ull*width*cells);s.after.resize(s.before.size());
        checked(cudaMemcpy(s.before.data(),reference,s.before.size()*sizeof(IntegrationReal),cudaMemcpyDeviceToHost));
    }
    return 1;
}
extern "C" int genesis_maint_state_consume(IntegrationReal *reference,IntegrationReal *first,IntegrationReal *second,int width,int cells,
    int short_copy,int work_copy,int long_copy) {
    auto &s=state;
    if(!s.pending)return 0;
    require(reference==s.reference&&width==s.width&&cells==s.cells,"reference storage changed during migration");
    require(!short_copy&&!work_copy&&!long_copy,"retired CPU state requested before GPU regeneration");
    if(genesis_maint_state_mode()==2) {
        const size_t bytes=3ull*width*cells*sizeof(IntegrationReal);
        checked(cudaMemcpy(s.after.data(),reference,bytes,cudaMemcpyDeviceToHost));
        require(std::memcmp(s.before.data(),s.after.data(),bytes)==0,"retained velocity reference changed");
        checked(cudaMemsetAsync(first,0xff,bytes,0));
        checked(cudaMemsetAsync(second,0xff,bytes,0));
        std::printf("GENESIS_MAINT_STATE: exact retained reference pass\n");
        std::printf("GENESIS_MAINT_STATE: GPU work poisoned\n");
    }
    s.pending=false;
    return 1;
}
