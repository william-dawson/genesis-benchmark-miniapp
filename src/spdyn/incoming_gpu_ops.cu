#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "migration_views.h"
#include "migration_coordinates.h"
using genesis::WorkingReal;
using genesis::IntegrationReal;
struct Records {
    const genesis::IntegrationReal *move_real,*stay_real;
    const int *move_int,*stay_int,*move_count,*stay_count;
    int move_real_rows,move_real_capacity,stay_real_rows,stay_real_capacity;
    int move_int_rows,move_int_capacity,stay_int_rows,stay_int_capacity;
};
extern "C" int genesis_state_cycle_publish(ParticleView,ParticleView,int);
extern "C" int genesis_pack_resident_boundary(ParticleView,ParticleView,const IntegrationReal*,IntegrationReal*,WorkingReal*);
extern "C" void genesis_water_payload_view(const IntegrationReal**,const IntegrationReal**,const IntegrationReal**,const int**,const int**,int*,int*);
extern "C" void genesis_hgr_payload_view(const genesis::IntegrationReal**,const int**,const int**,int*,int*,int*);
struct PackedView {const genesis::IntegrationReal *real=nullptr;const int *integer=nullptr,*ranges=nullptr;int sites=0;};
namespace {
void require(bool ok,const char *message) {
    if(!ok){std::fprintf(stderr,"GENESIS incoming GPU: %s\n",message);std::abort();}
}
void checked(cudaError_t code){require(code==cudaSuccess,cudaGetErrorString(code));}
int *error_device=nullptr;
void begin() {
    static const bool supported=[] {
        int pageable=0,ats=0;
        checked(cudaDeviceGetAttribute(&pageable,cudaDevAttrPageableMemoryAccess,0));
        checked(cudaDeviceGetAttribute(&ats,cudaDevAttrPageableMemoryAccessUsesHostPageTables,0));
        require(pageable&&ats,"coherent host page-table access is required");return true;
    }();
    (void)supported;
    if(!error_device)checked(cudaMalloc(&error_device,sizeof(int)));
    checked(cudaMemsetAsync(error_device,0,sizeof(int),0));
}
void finish() {
    checked(cudaGetLastError());
    int error=0;checked(cudaMemcpyAsync(&error,error_device,sizeof(int),cudaMemcpyDeviceToHost,0));
    checked(cudaStreamSynchronize(0));require(!error,"invalid incoming payload index or capacity");
}

int storage() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_INCOMING_STORAGE");if(!text)return 0;
        require(!std::strcmp(text,"0")||!std::strcmp(text,"1"),"storage must be 0 (coherent) or 1 (HBM)");
        return int(text[0]=='1');
    }();return mode;
}
bool oracle() {
    const char *value=std::getenv("GENESIS_INCOMING_GPU"),*legacy=std::getenv("GENESIS_MIGRATION_VERIFY");
    return (value&&value[0]=='2')||(legacy&&legacy[0]=='1');
}
struct CachedDomain {
    bool valid=false;
    ParticleView host,device;
    IntegrationReal *inverse=nullptr;
    WorkingReal *translation=nullptr,*water_charge=nullptr,*water_mass=nullptr;
    int *water_classes=nullptr;
} cache;
template<class T> T *clone(const T *host,size_t count) {
    T *device=nullptr;checked(cudaMalloc(&device,count*sizeof(T)));
    checked(cudaMemcpyAsync(device,host,count*sizeof(T),cudaMemcpyHostToDevice,0));return device;
}
void initialize(const ParticleView &p) {
    if(cache.valid) {require(cache.host.x==p.x&&cache.host.width==p.width&&cache.host.cells==p.cells,"cached domain shape changed");return;}
    require(p.noh&&p.groups&&p.hlocal,"cached domain needs constraint metadata");
    const IntegrationReal *x=nullptr,*v=nullptr,*packed=nullptr;const int *integer=nullptr,*ranges=nullptr;
    int width=0,cells=0;genesis_water_payload_view(&x,&v,&packed,&integer,&ranges,&width,&cells);
    require(x&&v&&width==p.width&&cells==p.cells,"resident MD state is unavailable");
    cache.host=p;cache.device=p;auto &d=cache.device;
    // Reuse the MD state after its outgoing consumers have completed packing.
    d.x=const_cast<IntegrationReal*>(x);d.v=const_cast<IntegrationReal*>(v);
    const size_t atoms=size_t(p.width)*p.cells;
    d.q=clone(p.q,atoms);d.m=clone(p.m,atoms);d.classes=clone(p.classes,atoms);
    d.ids=clone(p.ids,atoms);d.solute=clone(p.solute,atoms);d.reverse=clone(p.reverse,2ull*p.reverse_capacity);
    d.ns=clone(p.ns,p.cells);d.nw=clone(p.nw,p.cells);d.na=clone(p.na,p.cells);d.noh=clone(p.noh,p.cells);
    d.groups=clone(p.groups,size_t(p.group_height)*p.group_capacity*p.connect*p.cells);
    d.distance=clone(p.distance,size_t(p.dist_height)*p.dist_capacity*p.connect*p.cells);
    d.hlocal=clone(p.hlocal,size_t(p.connect)*p.cells);
    d.water=clone(p.water,size_t(p.water_height)*p.water_capacity*p.cells);
    checked(cudaMalloc(&cache.inverse,atoms*sizeof(IntegrationReal)));checked(cudaMalloc(&cache.translation,3*atoms*sizeof(WorkingReal)));
    cache.valid=true;
}
__global__ void publish_domain(ParticleView d,ParticleView host,int stage) {
    const int cell=blockIdx.x,ns=d.ns[cell],nw=d.nw[cell],noh=d.noh[cell];
    const int count=stage==0?noh:(stage==1?ns:d.na[cell]);
    for(int slot=threadIdx.x;slot<count;slot+=blockDim.x) {
        const int index=slot+d.width*cell;
        for(int axis=0;axis<3;++axis){host.x[3*index+axis]=d.x[3*index+axis];host.v[3*index+axis]=d.v[3*index+axis];}
        host.q[index]=d.q[index];host.m[index]=d.m[index];host.classes[index]=d.classes[index];host.ids[index]=d.ids[index];
        if(slot<(stage==0?noh:ns)) {
            const int id=d.solute[index];host.solute[index]=id;
            if(id>0&&id<=d.reverse_capacity){host.reverse[2*(id-1)]=d.reverse[2*(id-1)];host.reverse[2*(id-1)+1]=d.reverse[2*(id-1)+1];}
        }
    }
    if(!threadIdx.x) {host.noh[cell]=noh;if(stage>0)host.ns[cell]=ns;if(stage>1){host.nw[cell]=nw;host.na[cell]=d.na[cell];}}
    if(stage>0)for(int group=0;group<d.connect;++group) {
        const int count=d.hlocal[group+d.connect*cell],sites=group+2;
        if(!threadIdx.x)host.hlocal[group+d.connect*cell]=count;
        for(int i=threadIdx.x;i<sites*count;i+=blockDim.x) {
            const int entry=i/sites,site=i%sites;
            const size_t index=site+size_t(d.group_height)*(entry+d.group_capacity*(group+d.connect*cell));
            host.groups[index]=d.groups[index];
            const size_t real=site+size_t(d.dist_height)*(entry+d.dist_capacity*(group+d.connect*cell));
            host.distance[real]=d.distance[real];
        }
    }
    if(stage>1)for(int i=threadIdx.x;i<d.water_height*nw;i+=blockDim.x)
        host.water[i+size_t(d.water_height)*d.water_capacity*cell]=d.water[i+size_t(d.water_height)*d.water_capacity*cell];
}
void publish(int stage) {
    if(genesis_state_cycle_publish(cache.device,cache.host,stage))return;
    publish_domain<<<cache.device.cells,128>>>(cache.device,cache.host,stage);checked(cudaGetLastError());
}
__global__ void publish_boundary(ParticleView d,ParticleView host,IntegrationReal *inverse,WorkingReal *translation,IntegrationReal *host_inverse,WorkingReal *host_translation) {
    const int cell=blockIdx.x;
    for(int slot=threadIdx.x;slot<d.na[cell];slot+=blockDim.x) {
        const int index=slot+d.width*cell;host_inverse[index]=inverse[index];
        for(int axis=0;axis<3;++axis)host_translation[3*index+axis]=translation[3*index+axis];
        if(slot<d.ns[cell]) {
            const int id=d.solute[index];
            if(id>0&&id<=d.reverse_capacity){host.reverse[2*(id-1)]=d.reverse[2*(id-1)];host.reverse[2*(id-1)+1]=d.reverse[2*(id-1)+1];}
        }
    }
    if(!threadIdx.x)host.ns[cell]=d.ns[cell];
}
__device__ void bad(int *error){atomicExch(error,1);}
__device__ bool reverse(ParticleView p,int cell,int slot,int id,int *error) {
    if(id<1||id>p.reverse_capacity){bad(error);return false;}
    p.reverse[2*(id-1)]=static_cast<std::int16_t>(cell+1);
    p.reverse[2*(id-1)+1]=static_cast<std::int16_t>(slot+1);return true;
}
__device__ void atom_from_record(ParticleView p,int cell,int slot,const genesis::IntegrationReal *real,const int *integer,int *error) {
    const int index=slot+p.width*cell;
    for(int axis=0;axis<3;++axis){p.x[3*index+axis]=real[axis];p.v[3*index+axis]=real[axis+3];}
    p.q[index]=real[6];p.m[index]=real[7];p.classes[index]=integer[0];p.ids[index]=integer[1];p.solute[index]=integer[2];
    reverse(p,cell,slot,integer[2],error);
}
__device__ void atom_from_slot(ParticleView p,int cell,int destination,int source,int *error) {
    const int target=destination+p.width*cell,index=source+p.width*cell;
    for(int axis=0;axis<3;++axis){p.x[3*target+axis]=p.x[3*index+axis];p.v[3*target+axis]=p.v[3*index+axis];}
    p.q[target]=p.q[index];p.m[target]=p.m[index];p.classes[target]=p.classes[index];
    p.ids[target]=p.ids[index];p.solute[target]=p.solute[index];reverse(p,cell,destination,p.solute[target],error);
}
__global__ void incoming_atoms(ParticleView p,AtomRecords r,int *error) {
    const int cell=blockIdx.x*blockDim.x+threadIdx.x;
    if(cell>=p.cells)return;
    const int old=p.noh[cell],add=r.adds[cell],exits=r.exits[cell];
    if(old<0||old>p.width||add<0||exits<0||exits>old||old+add-exits>p.width||
       add>r.real_capacity||add>r.int_capacity||exits>r.index_capacity){bad(error);return;}
    for(int k=0;k<add;++k) {
        const int slot=k<exits?r.indices[k+r.index_capacity*cell]-1:k+old-exits;
        if(slot<0||slot>=p.width){bad(error);return;}
        atom_from_record(p,cell,slot,r.reals+r.real_rows*(k+r.real_capacity*cell),
            r.integers+r.int_rows*(k+r.int_capacity*cell),error);
    }
    if(add<exits) {
        int removed=0,source=old-1,next=add;
        while(removed<exits-add) {
            bool insert=true;
            for(int k=next;k<exits;++k)if(source==r.indices[k+r.index_capacity*cell]-1){insert=false;++removed;break;}
            if(insert) {
                const int slot=r.indices[next+r.index_capacity*cell]-1;
                if(source<0||slot<0||slot>=p.width){bad(error);return;}
                atom_from_slot(p,cell,slot,source,error);++removed;++next;
            }
            --source;
        }
    }
    p.noh[cell]=old+add-exits;
}
template<bool packed> __global__ void incoming_waters(ParticleView p,Records r,int sites,const WorkingReal *charge,const WorkingReal *mass,
    const int *classes,int *error,PackedView source_view) {
    const int cell=blockIdx.x;
    const int first_stay=packed?source_view.ranges[2*cell]:0,first_move=packed?source_view.ranges[2*(p.cells+cell)]:0;
    const int stays=packed?(first_stay<0?0:source_view.ranges[2*cell+1]-first_stay):r.stay_count[cell];
    const int moves=packed?(first_move<0?0:source_view.ranges[2*(p.cells+cell)+1]-first_move):r.move_count[cell];
    const int total=stays+moves;
    if(stays<0||moves<0||total>p.water_capacity||p.ns[cell]+sites*total>p.width||
        stays>r.stay_real_capacity||stays>r.stay_int_capacity||moves>r.move_real_capacity||moves>r.move_int_capacity) {
        if(!threadIdx.x)bad(error);return;
    }
    for(int i=threadIdx.x;i<sites*total;i+=blockDim.x) {
        const int entry=i/sites,site=i%sites,slot=p.ns[cell]+i,index=slot+p.width*cell;
        const bool stay=entry<stays;const int source=stay?entry:entry-stays;
        const int rr=stay?r.stay_real_rows:r.move_real_rows,rc=stay?r.stay_real_capacity:r.move_real_capacity;
        const int ir=stay?r.stay_int_rows:r.move_int_rows,ic=stay?r.stay_int_capacity:r.move_int_capacity;
        const auto *values=packed?source_view.real+18ull*((stay?first_stay:first_move)+source)+6*site:(stay?r.stay_real:r.move_real)+6*site+rr*(source+rc*cell);
        const int *ids=packed?source_view.integer+3ull*((stay?first_stay:first_move)+source)+site:(stay?r.stay_int:r.move_int)+site+ir*(source+ic*cell);
        for(int axis=0;axis<3;++axis){p.x[3*index+axis]=values[axis];p.v[3*index+axis]=values[axis+3];}
        p.ids[index]=*ids;p.q[index]=charge[site];p.m[index]=mass[site];p.classes[index]=classes[site];
        p.water[site+p.water_height*(entry+p.water_capacity*cell)]=slot+1;
    }
    if(!threadIdx.x){p.nw[cell]=total;p.na[cell]=p.ns[cell]+sites*total;}
}
template<bool packed> __global__ void incoming_groups(ParticleView p,Records r,int *error,PackedView source_view) {
    const int cell=blockIdx.x;
    int offset=p.noh[cell];
    for(int group=0;group<p.connect;++group) {
        const int stay_bucket=group+p.connect*cell,move_bucket=group+p.connect*(p.cells+cell);
        const int first_stay=packed?source_view.ranges[2*stay_bucket]:0,first_move=packed?source_view.ranges[2*move_bucket]:0;
        const int stays=packed?(first_stay<0?0:source_view.ranges[2*stay_bucket+1]-first_stay):r.stay_count[stay_bucket];
        const int moves=packed?(first_move<0?0:source_view.ranges[2*move_bucket+1]-first_move):r.move_count[stay_bucket];
        const int total=stays+moves,sites=group+2;
        if(stays<0||moves<0||total>p.group_capacity||total>p.dist_capacity||offset+sites*total>p.width||
           stays>r.stay_real_capacity||stays>r.stay_int_capacity||moves>r.move_real_capacity||moves>r.move_int_capacity) {
            if(!threadIdx.x)bad(error);return;
        }
        for(int i=threadIdx.x;i<sites*total;i+=blockDim.x) {
            const int entry=i/sites,site=i%sites,slot=offset+i,index=slot+p.width*cell;
            const bool stay=entry<stays;const int source=stay?entry:entry-stays;
            const int rr=stay?r.stay_real_rows:r.move_real_rows,rc=stay?r.stay_real_capacity:r.move_real_capacity;
            const int ir=stay?r.stay_int_rows:r.move_int_rows,ic=stay?r.stay_int_capacity:r.move_int_capacity;
            const auto *values=packed?source_view.real+9ull*source_view.sites*((stay?first_stay:first_move)+source)+9*site:(stay?r.stay_real:r.move_real)+9*site+rr*(source+rc*(group+p.connect*cell));
            const int *ids=packed?source_view.integer+3ull*source_view.sites*((stay?first_stay:first_move)+source)+3*site:(stay?r.stay_int:r.move_int)+3*site+ir*(source+ic*(group+p.connect*cell));
            for(int axis=0;axis<3;++axis){p.x[3*index+axis]=values[axis];p.v[3*index+axis]=values[axis+3];}
            p.q[index]=values[6];p.m[index]=values[7];p.classes[index]=ids[0];p.ids[index]=ids[1];p.solute[index]=ids[2];
            p.groups[site+p.group_height*(entry+p.group_capacity*(group+p.connect*cell))]=slot+1;
            p.distance[site+p.dist_height*(entry+p.dist_capacity*(group+p.connect*cell))]=values[8];
            reverse(p,cell,slot,ids[2],error);
        }
        if(!threadIdx.x)p.hlocal[group+p.connect*cell]=total;
        offset+=sites*total;
    }
    if(!threadIdx.x)p.ns[cell]=offset;
}
struct Geometry {double box[3],origin[3];int sites;double epsilon;};
__device__ double periodic_translation(double x,double box,double origin) {
#if defined(_SINGLE)
    return genesis::single_translation(x,origin,box);
#else
    const double shift=__dsub_rn(x,origin);
    const double move=__dsub_rn(__dmul_rn(box,0.5),__dmul_rn(box,round(__ddiv_rn(shift,box))));
    return __dsub_rn(move,origin);
#endif
}
__device__ void translate(ParticleView p,int cell,int source,int destination,WorkingReal *translation,Geometry g,int *error,bool solute) {
    for(int axis=0;axis<3;++axis)
        translation[3*(destination+p.width*cell)+axis]=periodic_translation(p.x[3*(source+p.width*cell)+axis],g.box[axis],g.origin[axis]);
    if(solute)reverse(p,cell,destination,p.solute[destination+p.width*cell],error);
}
__global__ void cell_bookkeeping(ParticleView p,IntegrationReal *inverse,WorkingReal *translation,Geometry g,int *error) {
    const int cell=blockIdx.x,count=p.na[cell],noh=p.noh[cell];
    if(count<0||count>p.width||noh<0||noh>count){if(!threadIdx.x)bad(error);return;}
    for(int slot=threadIdx.x;slot<count;slot+=blockDim.x) {
        const int index=slot+p.width*cell;
#if defined(_SINGLE)
        inverse[index]=p.m[index]>g.epsilon?__fdiv_rn(1.0f,p.m[index]):0.0f;
#else
        inverse[index]=p.m[index]>g.epsilon?__ddiv_rn(1.0,p.m[index]):0.0;
#endif
    }
    for(int slot=threadIdx.x;slot<noh;slot+=blockDim.x)translate(p,cell,slot,slot,translation,g,error,true);
    int offset=noh;
    for(int group=0;group<p.connect;++group) {
        const int groups=p.hlocal[group+p.connect*cell],sites=group+2;
        if(groups<0||groups>p.group_capacity||offset+groups*sites>count){if(!threadIdx.x)bad(error);return;}
        for(int entry=threadIdx.x;entry<groups;entry+=blockDim.x) {
            const int source=p.groups[p.group_height*(entry+p.group_capacity*(group+p.connect*cell))]-1;
            if(source<0||source>=count){bad(error);continue;}
            double move[3];
            for(int axis=0;axis<3;++axis)move[axis]=periodic_translation(p.x[3*(source+p.width*cell)+axis],g.box[axis],g.origin[axis]);
            for(int site=0;site<sites;++site) {
                const int slot=offset+sites*entry+site;
                for(int axis=0;axis<3;++axis)translation[3*(slot+p.width*cell)+axis]=move[axis];
                reverse(p,cell,slot,p.solute[slot+p.width*cell],error);
            }
        }
        offset+=sites*groups;
    }
    if(!threadIdx.x)p.ns[cell]=offset;
    if((count-offset)%g.sites){if(!threadIdx.x)bad(error);return;}
    const int waters=(count-offset)/g.sites;
    for(int entry=threadIdx.x;entry<waters;entry+=blockDim.x) {
        const int first=offset+g.sites*entry;
        double move[3];
        for(int axis=0;axis<3;++axis)move[axis]=periodic_translation(p.x[3*(first+p.width*cell)+axis],g.box[axis],g.origin[axis]);
        for(int site=0;site<g.sites;++site)for(int axis=0;axis<3;++axis)
            translation[3*(first+site+p.width*cell)+axis]=move[axis];
    }
}

}
extern "C" int genesis_incoming_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_INCOMING_GPU");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");return int(value);
    }();return mode;
}
extern "C" int genesis_atom_packets(AtomRecords*);
extern "C" void genesis_incoming_atoms(const ParticleView *view,const AtomRecords *records) {
    AtomRecords packets=*records;genesis_atom_packets(&packets);records=&packets;
    begin();
    if(storage()) {
        initialize(*view);incoming_atoms<<<(view->cells+127)/128,128>>>(cache.device,*records,error_device);
        if(oracle())publish(0);
    } else incoming_atoms<<<(view->cells+127)/128,128>>>(*view,*records,error_device);
    finish();
}
extern "C" void genesis_incoming_waters(const ParticleView *view,const Records *records,int sites,
    const WorkingReal *charge,const WorkingReal *mass,const int *classes) {
    begin();PackedView packed;
    if(storage()) {
        require(cache.valid&&sites==3,"cached water requires resident TIP3 state");
        const IntegrationReal *x=nullptr,*v=nullptr;int width=0,cells=0;
        genesis_water_payload_view(&x,&v,&packed.real,&packed.integer,&packed.ranges,&width,&cells);
        require(packed.real&&packed.integer&&packed.ranges&&width==view->width&&cells==view->cells,"missing water payload");
        if(!cache.water_charge){cache.water_charge=clone(charge,sites);cache.water_mass=clone(mass,sites);cache.water_classes=clone(classes,sites);}
        incoming_waters<true><<<view->cells,128>>>(cache.device,*records,sites,cache.water_charge,cache.water_mass,cache.water_classes,error_device,packed);
        publish(2);
    } else incoming_waters<false><<<view->cells,128>>>(*view,*records,sites,charge,mass,classes,error_device,packed);
    finish();
}
extern "C" void genesis_incoming_groups(const ParticleView *view,const Records *records) {
    begin();PackedView packed;
    if(storage()) {
        require(cache.valid,"cached group state is missing");int cells=0,connect=0;
        genesis_hgr_payload_view(&packed.real,&packed.integer,&packed.ranges,&packed.sites,&cells,&connect);
        require(packed.real&&packed.integer&&packed.ranges&&cells==view->cells&&connect==view->connect,"missing group payload");
        incoming_groups<true><<<view->cells,128>>>(cache.device,*records,error_device,packed);
        if(oracle())publish(1);
    } else incoming_groups<false><<<view->cells,128>>>(*view,*records,error_device,packed);
    finish();
}
extern "C" void genesis_incoming_boundary(const ParticleView *view,IntegrationReal *inverse,WorkingReal *translation,
    const double *box,const double *origin,int sites,double epsilon) {
    require(sites>=1&&sites<=4,"unsupported water size");
    Geometry geometry;
    for(int axis=0;axis<3;++axis){require(box[axis]>0,"invalid box");geometry.box[axis]=box[axis];geometry.origin[axis]=origin[axis];}
    geometry.sites=sites;geometry.epsilon=epsilon;
    begin();
    if(storage()&&cache.valid) {
        cell_bookkeeping<<<view->cells,128>>>(cache.device,cache.inverse,cache.translation,geometry,error_device);
        if(!genesis_pack_resident_boundary(cache.device,cache.host,cache.inverse,inverse,translation))
          publish_boundary<<<view->cells,128>>>(cache.device,cache.host,cache.inverse,cache.translation,inverse,translation);
    } else cell_bookkeeping<<<view->cells,128>>>(*view,inverse,translation,geometry,error_device);
    finish();
}

extern "C" int genesis_incoming_cached_view(const ParticleView *host,ParticleView *device) {
 if(!storage())return 0;initialize(*host);*device=cache.device;return 1;
}
extern "C" int genesis_incoming_pair_view(ParticleView *view,const WorkingReal **translation) {
 if(!storage() || !cache.valid)return 0;*view=cache.device;*translation=cache.translation;return 1;
}
extern "C" int genesis_incoming_metadata_view(ParticleView *view,const IntegrationReal **inverse) {
 if(!storage() || !cache.valid)return 0;*view=cache.device;*inverse=cache.inverse;return 1;
}
