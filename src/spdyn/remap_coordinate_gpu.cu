#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include "precision_types.h"
using genesis::WorkingReal;
using genesis::IntegrationReal;

namespace {
bool pending=false;
unsigned consumers=0;
void require(bool condition,const char *message) {
    if(!condition) {std::fprintf(stderr,"GENESIS remap coordinates: %s\n",message);std::abort();}
}
void checked(cudaError_t code) {
    if(code!=cudaSuccess) {std::fprintf(stderr,"GENESIS remap coordinates: %s\n",cudaGetErrorString(code));std::abort();}
}
bool selected(const char *name) {
    const char *value=std::getenv(name);return value&&std::strcmp(value,"1")==0;
}
}
extern "C" int genesis_remap_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_REMAP_PREPARE");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");return int(value);
    }();return mode;
}
extern "C" int genesis_remap_allowed() {
    return genesis_remap_mode()&&selected("GENESIS_FORCE_BONDED");
}
extern "C" void genesis_remap_poison(WorkingReal *packed,double *raw,int atoms) {
    if(genesis_remap_mode()!=2)return;
    require(!pending,"unconsumed prior publication");
    checked(cudaMemsetAsync(packed,0xff,4ull*atoms*sizeof(WorkingReal),0));
    require(raw,"raw PME buffer missing");checked(cudaMemsetAsync(raw,0xff,4ull*atoms*sizeof(double),0));
    std::printf("GENESIS_REMAP: GPU coordinate packs poisoned\n");
}
extern "C" void genesis_remap_publish(const IntegrationReal *coordinates,const IntegrationReal *velocity,const WorkingReal *translation,
    const WorkingReal *packed,const double *raw,const IntegrationReal *host_coordinates,const IntegrationReal *host_velocity,
    const WorkingReal *host_translation,const WorkingReal *charge,const int *counts,const int *starts,int width,int cells,int atoms) {
    require(!pending&&genesis_remap_mode(),"invalid remap publication");
    consumers=0;pending=true;
    if(genesis_remap_mode()!=2)return;
    std::vector<IntegrationReal> x(3ull*width*cells),v(x.size());
    std::vector<double> q(4ull*atoms);
    std::vector<WorkingReal> t(x.size()),p(4ull*atoms);
    checked(cudaMemcpy(x.data(),coordinates,x.size()*sizeof(IntegrationReal),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(v.data(),velocity,v.size()*sizeof(IntegrationReal),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(t.data(),translation,t.size()*sizeof(WorkingReal),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(p.data(),packed,p.size()*sizeof(WorkingReal),cudaMemcpyDeviceToHost));
    checked(cudaMemcpy(q.data(),raw,q.size()*sizeof(double),cudaMemcpyDeviceToHost));
    int total=0;
    for(int cell=0;cell<cells;++cell) {
        require(counts[cell]>=0&&counts[cell]<=width&&starts[cell]==total,"invalid remap offsets");
        for(int atom=0;atom<counts[cell];++atom) {
            const size_t local=size_t(cell)*width+atom,global=starts[cell]+atom;
            for(int axis=0;axis<3;++axis) {
                const size_t i=3*local+axis;
                require(std::isfinite(host_coordinates[i])&&std::isfinite(host_velocity[i])&&std::isfinite(host_translation[i]),"nonfinite remap source");
                const IntegrationReal expected_x=host_coordinates[i],expected_v=host_velocity[i];
                require(std::memcmp(x.data()+i,&expected_x,sizeof(IntegrationReal))==0,"raw coordinate mismatch");
                require(std::memcmp(v.data()+i,&expected_v,sizeof(IntegrationReal))==0,"velocity mismatch");
                require(std::memcmp(t.data()+i,host_translation+i,sizeof(WorkingReal))==0,"unit scaling changed translation");
                const WorkingReal expected=host_coordinates[i]+host_translation[i];
                require(std::memcmp(p.data()+size_t(axis)*atoms+global,&expected,sizeof(WorkingReal))==0,"periodic coordinate mismatch");
                const double expected_raw=expected_x;
                require(std::memcmp(q.data()+3*global+axis,&expected_raw,sizeof(double))==0,"PME raw coordinate mismatch");
            }
            require(std::memcmp(p.data()+3ull*atoms+global,charge+local,sizeof(WorkingReal))==0,"periodic charge mismatch");
            const double expected_charge=charge[local];
            require(std::memcmp(q.data()+3ull*atoms+global,&expected_charge,sizeof(double))==0,"PME charge mismatch");
        }
        total+=counts[cell];
    }
    require(total==atoms,"atom count mismatch");
    std::printf("GENESIS_REMAP: exact coordinates velocity translations and packs pass\n");
}
extern "C" void genesis_remap_record(int consumer) {
    if(pending) {require(consumer>=0&&consumer<=3,"unknown consumer");consumers|=1u<<consumer;}
}
extern "C" void genesis_remap_consume() {
    if(!pending)return;
    require((consumers&14u)==14u,"bonded raw/periodic or PME input bypassed remap publication");
    pending=false;
    if(genesis_remap_mode()==2)std::printf("GENESIS_REMAP: all required GPU consumers pass\n");
}
