#include <cuda_runtime.h>
#include "../config.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

extern "C" void genesis_translation_require_host();
extern "C" void genesis_translation_invalidate();
extern "C" void genesis_translation_configure(const double*,const double*,const int*,int,int);
extern "C" void genesis_translation_configure_device(const double*,const double*,const int*,int,int);
extern "C" int genesis_pack_resident_available();
extern "C" int genesis_pack_resident_copy(const double*,double*,double*,int*,int*,const double*,const double*,const int*,const int*,int,int);
extern "C" int genesis_state_cycle_retired();
extern "C" void genesis_state_cycle_host_refreshed();

extern "C" int genesis_nonbond_virial_host_required();
extern "C" int genesis_remap_mode();
extern "C" void genesis_remap_poison(double*,double*,int);
extern "C" void genesis_remap_record(int);
extern "C" void genesis_remap_publish(const double*,const double*,const double*,const double*,const double*,
    const double*,const double*,const double*,const double*,const int*,const int*,int,int,int);

extern "C" int genesis_pair_prepare_coord_copy(void*,int,int,std::uintptr_t,int);
extern "C" void genesis_pair_prepare_expire();
namespace {
void checked(cudaError_t status) {
    if(status!=cudaSuccess) {
        std::fprintf(stderr,"GENESIS coordinate packing: %s\n",cudaGetErrorString(status));
        std::abort();
    }
}
int setting(const char *name,int maximum=1) {
    const char *value=std::getenv(name);
    if(!value)return 0;
    char *end=nullptr;
    const long selected=std::strtol(value,&end,10);
    if(end==value||*end||selected<0||selected>maximum) {
        std::fprintf(stderr,"Invalid %s\n",name);std::abort();
    }
    return int(selected);
}
struct State {
    const double *coordinates=nullptr;
    double *translation=nullptr,*charge=nullptr,*packed=nullptr,*raw=nullptr;
    int *counts=nullptr,*starts=nullptr;
    int width=0,cells=0,atoms=0;
    bool active=false,short_step=false,host_current=true,copied[3]{};
    cudaEvent_t ready=nullptr,consumed[3]{};
} state;

__global__ void pack_periodic(const double *coordinates,double *translation,const double *charge,
    const int *counts,const int *starts,double *packed,double *raw,int width,int cells,int atoms,double3 scale) {
    const int index=blockIdx.x*blockDim.x+threadIdx.x;
    const int cell=index/width,atom=index%width;
    if(cell>=cells||atom>=counts[cell])return;
    const int global=starts[cell]+atom;
    const size_t base=3ull*index;
    const double tx=__dmul_rn(translation[base],scale.x);
    const double ty=__dmul_rn(translation[base+1],scale.y);
    const double tz=__dmul_rn(translation[base+2],scale.z);
    translation[base]=tx;translation[base+1]=ty;translation[base+2]=tz;
    packed[global]=__dadd_rn(coordinates[base],tx);
    packed[atoms+global]=__dadd_rn(coordinates[base+1],ty);
    packed[2ull*atoms+global]=__dadd_rn(coordinates[base+2],tz);
    packed[3ull*atoms+global]=charge[index];
    if(raw) {
        raw[3ull*global]=coordinates[base];
        raw[3ull*global+1]=coordinates[base+1];
        raw[3ull*global+2]=coordinates[base+2];
        raw[3ull*atoms+global]=charge[index];
    }
}
}

extern "C" int genesis_coord_enabled() {
    static const int enabled=setting("GENESIS_COORD_GPU");
#if defined(_SINGLE) || defined(_MIXED)
    if(enabled)std::abort();
#endif
    return enabled;
}
extern "C" int genesis_coord_active() {return state.active;}
extern "C" int genesis_coord_resident() {
    static const int enabled=setting("GENESIS_COORD_RESIDENT");return enabled;
}
extern "C" int genesis_coord_verify() {
    static const int enabled=setting("GENESIS_COORD_VERIFY",2);return enabled;
}
extern "C" int genesis_coord_pme() {
    static const int enabled=setting("GENESIS_COORD_PME");
    if(enabled&&(!genesis_coord_enabled()||!setting("GENESIS_PME_GPU")||!setting("GENESIS_PME_PACKED"))) {
        std::fprintf(stderr,"Resident PME input requires GPU coordinates and packed GPU PME\n");std::abort();
    }
    return enabled;
}
extern "C" void genesis_coord_host_refreshed() {state.host_current=true;genesis_state_cycle_host_refreshed();}
extern "C" void genesis_coord_require_host() {
    genesis_translation_require_host();
    if(state.active&&!state.host_current) {
        std::fprintf(stderr,"A CPU packed-coordinate consumer has stale coordinates\n");std::abort();
    }
}
extern "C" void genesis_coord_invalidate() {genesis_translation_invalidate();state.active=false;}

extern "C" void genesis_coord_configure(const double *coordinates,const double *translation,
    const double *charge,const int *counts,const int *starts,int width,int cells,int atoms) {
    if(!genesis_coord_enabled())return;
    auto &s=state;
    if(width<=0||cells<=0||atoms<=0)std::abort();
    if(s.width==0) {
        s.width=width;s.cells=cells;s.atoms=atoms;
        checked(cudaMalloc(&s.translation,3ull*width*cells*sizeof(double)));
        checked(cudaMalloc(&s.charge,size_t(width)*cells*sizeof(double)));
        checked(cudaMalloc(&s.packed,4ull*atoms*sizeof(double)));
        checked(cudaMalloc(&s.raw,4ull*atoms*sizeof(double)));
        checked(cudaMalloc(&s.counts,size_t(cells)*sizeof(int)));
        checked(cudaMalloc(&s.starts,size_t(cells)*sizeof(int)));
        checked(cudaEventCreateWithFlags(&s.ready,cudaEventDisableTiming));
        for(auto &event:s.consumed)checked(cudaEventCreateWithFlags(&event,cudaEventDisableTiming));
    }
    if(width!=s.width||cells!=s.cells||atoms!=s.atoms)std::abort();
    int next=0;
    for(int cell=0;cell<cells;++cell) {
        if(counts[cell]<0||counts[cell]>width||starts[cell]!=next)std::abort();
        next+=counts[cell];
    }
    if(next!=atoms)std::abort();
    s.coordinates=coordinates;s.active=false;
    if(genesis_pack_resident_copy(coordinates,s.translation,s.charge,s.counts,s.starts,translation,charge,counts,starts,width,cells)) {
        genesis_translation_configure_device(s.translation,translation,counts,width,cells);return;
    }
    genesis_translation_configure(s.translation,translation,counts,width,cells);
    checked(cudaMemcpy(s.translation,translation,3ull*width*cells*sizeof(double),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.charge,charge,size_t(width)*cells*sizeof(double),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.counts,counts,size_t(cells)*sizeof(int),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(s.starts,starts,size_t(cells)*sizeof(int),cudaMemcpyHostToDevice));
}

extern "C" void genesis_coord_pack(double x,double y,double z,int short_step) {
    genesis_pair_prepare_expire();
    if(!genesis_coord_enabled())return;
    auto &s=state;
    if(!s.coordinates)std::abort();
    for(int consumer=0;consumer<3;++consumer)if(s.copied[consumer]) {
        checked(cudaStreamWaitEvent(0,s.consumed[consumer],0));s.copied[consumer]=false;
    }
    pack_periodic<<<(s.width*s.cells+255)/256,256>>>(s.coordinates,s.translation,s.charge,s.counts,s.starts,
        s.packed,genesis_coord_pme()?s.raw:nullptr,s.width,s.cells,s.atoms,make_double3(x,y,z));
    checked(cudaGetLastError());checked(cudaEventRecord(s.ready,0));
    s.active=true;s.short_step=short_step!=0;s.host_current=false;
}

extern "C" int genesis_coord_copy(void *destination,int atoms,int components,std::uintptr_t stream,int consumer) {
    auto &s=state;
    if(!s.active)return genesis_pair_prepare_coord_copy(destination,atoms,components,stream,consumer);
    if(atoms!=s.atoms||components<3||components>4||consumer<0||consumer>2)std::abort();
    if(consumer==2&&!genesis_coord_pme())return 0;
    auto target=reinterpret_cast<cudaStream_t>(stream);
    checked(cudaStreamWaitEvent(target,s.ready,0));
    checked(cudaMemcpyAsync(destination,consumer==2?s.raw:s.packed,size_t(components)*atoms*sizeof(double),cudaMemcpyDeviceToDevice,target));
    checked(cudaEventRecord(s.consumed[consumer],target));s.copied[consumer]=true;
    genesis_remap_record(consumer);
    return 1;
}

extern "C" int genesis_coord_cpu_required(int virial,int cpu_calc) {
    if(!state.active)return 1;
    const bool required=genesis_coord_verify()==1||setting("GENESIS_FORCE_VERIFY")||(!state.short_step&&(cpu_calc||(virial&&genesis_nonbond_virial_host_required())));
    if(required)genesis_coord_require_host();
    return required;
}

extern "C" void genesis_coord_check(const double *coordinates,const double *translation,const double *charge,
    const int *counts,const int *starts) {
    const int mode=genesis_coord_verify();
    if(!mode||!state.active)return;
    auto &s=state;
    checked(cudaEventSynchronize(s.ready));
    std::vector<double> packed(4ull*s.atoms),snapshot;
    checked(cudaMemcpy(packed.data(),s.packed,packed.size()*sizeof(double),cudaMemcpyDeviceToHost));
    if(mode==2) {
        snapshot.resize(3ull*s.width*s.cells);
        checked(cudaMemcpy(snapshot.data(),s.coordinates,snapshot.size()*sizeof(double),cudaMemcpyDeviceToHost));
        coordinates=snapshot.data();
    }
    std::vector<double> raw;
    if(genesis_coord_pme()) {
        raw.resize(4ull*s.atoms);
        checked(cudaMemcpy(raw.data(),s.raw,raw.size()*sizeof(double),cudaMemcpyDeviceToHost));
    }
    double maximum=0;
    for(int cell=0;cell<s.cells;++cell)for(int atom=0;atom<counts[cell];++atom) {
        const size_t local=size_t(cell)*s.width+atom;
        const size_t global=starts[cell]+atom;
        for(int axis=0;axis<4;++axis) {
            const double expected=axis==3?charge[local]:coordinates[3*local+axis]+translation[3*local+axis];
            const double actual=packed[size_t(axis)*s.atoms+global];
            if(!std::isfinite(actual)||!std::isfinite(expected))std::abort();
            const double ratio=std::abs(actual-expected)/(axis==3?1e-15:1e-10+1e-12*std::abs(expected));
            maximum=std::fmax(maximum,ratio);
            if(!raw.empty()) {
                const double expected_raw=axis==3?charge[local]:coordinates[3*local+axis];
                const double actual_raw=axis==3?raw[3ull*s.atoms+global]:raw[3ull*global+axis];
                if(!std::isfinite(actual_raw)||actual_raw!=expected_raw)std::abort();
            }
        }
    }
    if(maximum>1) {std::fprintf(stderr,"GENESIS periodic packing error %.12e\n",maximum);std::abort();}
    std::printf("GENESIS_COORD_VERIFY: periodic coordinates pass %.12e\n",maximum);
}

extern "C" void genesis_coord_prepare_remap(const double *coordinates,const double *velocity,const double *gpu_velocity,
    const double *translation,const double *charge,const int *counts,const int *starts) {
    auto &s=state;
    if(!genesis_remap_mode()||s.active||!s.coordinates||!genesis_coord_pme())std::abort();
    if(genesis_remap_mode()==2) {
        for(int consumer=0;consumer<3;++consumer)if(s.copied[consumer])
            checked(cudaStreamWaitEvent(0,s.consumed[consumer],0));
        genesis_remap_poison(s.packed,s.raw,s.atoms);
    }
    genesis_coord_pack(1.0,1.0,1.0,0);
    s.host_current=!genesis_state_cycle_retired();
    genesis_remap_publish(s.coordinates,gpu_velocity,s.translation,s.packed,s.raw,coordinates,velocity,
        translation,charge,counts,starts,s.width,s.cells,s.atoms);
}
