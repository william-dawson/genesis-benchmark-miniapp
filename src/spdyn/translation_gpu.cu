#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {
struct State {
    const double *device=nullptr,*host=nullptr;
    int width=0,cells=0;
    bool dirty=false,rebuilding=false,rebuilt_device=false;
    std::vector<int> counts;
    std::vector<double> observed;
} state;
void checked(cudaError_t code) {
    if(code!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS translation: %s\n",cudaGetErrorString(code));std::abort();
    }
}
void require(bool condition,const char *message) {
    if(!condition) {std::fprintf(stderr,"GENESIS translation: %s\n",message);std::abort();}
}
}
extern "C" int genesis_translation_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_TRANSLATION_GPU");
        if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");
        return int(value);
    }();
    return mode;
}
extern "C" int genesis_translation_rebuild_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_TRANSLATION_REBUILD");
        if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"rebuild mode must be 0, 1 or 2");
        return int(value);
    }();
    return mode;
}
extern "C" int genesis_translation_rebuild_begin() {
    if(!genesis_translation_rebuild_mode())return 0;
    require(genesis_translation_mode()!=0,"rebuild retirement requires GPU translations");
    if(!state.device)return 0;
    require(!state.rebuilding,"nested translation rebuild");
    state.dirty=false;state.rebuilding=true;
    return 1;
}
extern "C" void genesis_translation_poisoned() {
    require(state.device&&!state.rebuilding,"invalid output mirror poison");
    state.dirty=true;
}
extern "C" int genesis_translation_rebuilding() {return state.rebuilding;}
extern "C" void genesis_translation_rebuild_end() {
    require(state.rebuilding&&!state.dirty,"invalid translation rebuild completion");
    state.rebuilding=false;
    state.dirty=state.rebuilt_device;state.rebuilt_device=false;
}
extern "C" void genesis_translation_rebuilt_device() {
    require(state.rebuilding,"device translation rebuild was not declared");state.rebuilt_device=true;
}
extern "C" void genesis_translation_invalidate() {
    if(genesis_translation_mode())require(!state.dirty,"invalidation needs refresh or a declared rebuild");
}
extern "C" void genesis_translation_require_host() {
    if(genesis_translation_mode())require(!state.dirty&&!state.rebuilding,"CPU consumer or reconfiguration has a stale mirror");
}
void configure(const double *device,const double *host,const int *counts,int width,int cells,bool resident) {
    if(!genesis_translation_mode())return;
    if(!resident)genesis_translation_require_host();
    else require(!state.rebuilding,"device reconfiguration during rebuild");
    require(device&&host&&width>0&&cells>0,"invalid configuration");
    auto &s=state;
    if(s.device)require(device==s.device&&host==s.host&&width==s.width&&cells==s.cells,"translation layout changed");
    s.device=device;s.host=host;s.width=width;s.cells=cells;
    s.counts.assign(counts,counts+cells);
    if(genesis_translation_mode()==2)s.observed.resize(3ull*width*cells);
}
extern "C" void genesis_translation_configure(const double *device,const double *host,const int *counts,int width,int cells) {
    configure(device,host,counts,width,cells,false);
}
extern "C" void genesis_translation_configure_device(const double *device,const double *host,const int *counts,int width,int cells) {
    configure(device,host,counts,width,cells,true);
}
extern "C" void genesis_translation_publish(const double *reference) {
    if(!genesis_translation_mode())return;
    auto &s=state;
    require(s.device&&!s.rebuilding,"publication before configuration or rebuild completion");
    s.dirty=true;
    if(genesis_translation_mode()!=2)return;
    checked(cudaMemcpy(s.observed.data(),s.device,s.observed.size()*sizeof(double),cudaMemcpyDeviceToHost));
    for(int cell=0;cell<s.cells;++cell)for(int atom=0;atom<s.counts[cell];++atom) {
        const size_t offset=3ull*(size_t(cell)*s.width+atom);
        for(int axis=0;axis<3;++axis) {
            require(std::isfinite(reference[offset+axis])&&std::isfinite(s.observed[offset+axis]),"nonfinite translation");
            require(std::memcmp(reference+offset+axis,s.observed.data()+offset+axis,sizeof(double))==0,"CPU/GPU scaling differs");
        }
    }
    std::printf("GENESIS_TRANSLATION_VERIFY: exact scaling pass\n");
}
extern "C" void genesis_translation_refresh(double *host) {
    if(!genesis_translation_mode())return;
    auto &s=state;
    if(s.rebuilding||!s.dirty)return;
    require(host==s.host,"refresh target changed");
    checked(cudaMemcpy(host,s.device,3ull*s.width*s.cells*sizeof(double),cudaMemcpyDeviceToHost));
    s.dirty=false;
    if(genesis_translation_mode()==2)std::printf("GENESIS_TRANSLATION_VERIFY: CPU mirror refreshed\n");
}
