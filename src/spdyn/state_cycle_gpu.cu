#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>
#include "migration_views.h"

extern "C" int genesis_incoming_metadata_view(ParticleView*,const double**);
extern "C" int genesis_incoming_pair_view(ParticleView*,const double**);
extern "C" int genesis_exclusion_prefix(const int**,const int**);
extern "C" void genesis_translation_rebuilt_device();
namespace {
void require(bool ok,const char *s){if(!ok){std::fprintf(stderr,"GENESIS state cycle: %s\n",s);std::abort();}}
void checked(cudaError_t code){require(code==cudaSuccess,cudaGetErrorString(code));}
int setting(const char *name){const char *s=std::getenv(name);if(!s)return 0;
 require(std::strlen(s)==1&&s[0]>='0'&&s[0]<='3',name);return s[0]-'0';}
int state_mode(){static int mode=setting("GENESIS_STATE_CYCLE");return mode;}
int pack_mode(){static int mode=setting("GENESIS_PACK_RESIDENT");return mode;}
bool production(){
 for(const char *s:{"GENESIS_FORCE_VERIFY","GENESIS_STEP_VERIFY","GENESIS_STATE_VERIFY","GENESIS_THERMO_VERIFY",
  "GENESIS_COORD_VERIFY","GENESIS_WATER_VERIFY","GENESIS_HGR_VERIFY","GENESIS_MIGRATION_VERIFY",
  "GENESIS_BONDED_MIGRATION_VERIFY","GENESIS_EXCLUSION_VERIFY","GENESIS_PAIR_VERIFY",
  "GENESIS_PAIR_SORT_VERIFY","GENESIS_PAIR_COMPACT_VERIFY"})if(setting(s))return false;
 for(const char *s:{"GENESIS_INCOMING_GPU","GENESIS_INCOMING_STORAGE","GENESIS_ATOM_OUTGOING_GPU",
  "GENESIS_MIGRATION_RESIDENT","GENESIS_PAIR_PREPARE_GPU","GENESIS_EXCLUSION_GPU",
  "GENESIS_FORCE_TOPOLOGY_GPU","GENESIS_TRANSLATION_GPU","GENESIS_TRANSLATION_REBUILD",
  "GENESIS_REMAP_PREPARE","GENESIS_MAINT_STATE","GENESIS_LONG_FORCE","GENESIS_BONDED_TOTALS",
  "GENESIS_NONBOND_VIRIAL","GENESIS_NONBOND_MIRROR","GENESIS_FINAL_REDUCTIONS",
  "GENESIS_BONDED_INCOMING_GPU","GENESIS_BONDED_RESIDENCY","GENESIS_IMPROPER_OUTGOING_GPU"})if(setting(s)!=1)return false;
 return true;
}
struct Cycle {bool pending=false,retired=false;const double *x=nullptr,*v=nullptr;int width=0,cells=0;} cycle;
void compare_vectors(const double *device,const double *host,const int *counts,int width,int cells,bool poison){
 std::vector<double> values(3ull*width*cells);checked(cudaMemcpy(values.data(),device,values.size()*8,cudaMemcpyDeviceToHost));
 for(int c=0;c<cells;++c){require(counts[c]>=0&&counts[c]<=width,"invalid active count");
 for(int a=0;a<counts[c];++a)for(int k=0;k<3;++k){size_t i=3ull*(a+width*c)+k;
  require(std::isfinite(values[i]),"nonfinite resident state");
  require(poison?std::isnan(host[i]):std::memcmp(values.data()+i,host+i,8)==0,"state differs or CPU mirror was rewritten");}}
}
__global__ void publish_particle_metadata(ParticleView d,ParticleView h,int stage){
 int c=blockIdx.x,ns=d.ns[c],nw=d.nw[c],noh=d.noh[c],n=stage==0?noh:(stage==1?ns:d.na[c]);
 for(int a=threadIdx.x;a<n;a+=blockDim.x){int i=a+d.width*c;
  h.q[i]=d.q[i];h.m[i]=d.m[i];h.classes[i]=d.classes[i];h.ids[i]=d.ids[i];
  if(a<(stage==0?noh:ns)){int id=d.solute[i];h.solute[i]=id;
   if(id>0&&id<=d.reverse_capacity){h.reverse[2*(id-1)]=d.reverse[2*(id-1)];h.reverse[2*(id-1)+1]=d.reverse[2*(id-1)+1];}}
 }
 if(!threadIdx.x){h.noh[c]=noh;if(stage>0)h.ns[c]=ns;if(stage>1){h.nw[c]=nw;h.na[c]=d.na[c];}}
 if(stage>0)for(int g=0;g<d.connect;++g){int count=d.hlocal[g+d.connect*c],sites=g+2;
  if(!threadIdx.x)h.hlocal[g+d.connect*c]=count;
  for(int i=threadIdx.x;i<sites*count;i+=blockDim.x){int j=i/sites,a=i%sites;
   size_t k=a+size_t(d.group_height)*(j+d.group_capacity*(g+d.connect*c));h.groups[k]=d.groups[k];
   k=a+size_t(d.dist_height)*(j+d.dist_capacity*(g+d.connect*c));h.distance[k]=d.distance[k];}}
 if(stage>1)for(int i=threadIdx.x;i<d.water_height*nw;i+=blockDim.x)
  h.water[i+size_t(d.water_height)*d.water_capacity*c]=d.water[i+size_t(d.water_height)*d.water_capacity*c];
}
__global__ void publish_boundary_metadata(ParticleView d,ParticleView h,const double *inverse,double *host_inverse){
 int c=blockIdx.x;for(int a=threadIdx.x;a<d.na[c];a+=blockDim.x){int i=a+d.width*c;host_inverse[i]=inverse[i];
  if(a<d.ns[c]){int id=d.solute[i];if(id>0&&id<=d.reverse_capacity){h.reverse[2*(id-1)]=d.reverse[2*(id-1)];h.reverse[2*(id-1)+1]=d.reverse[2*(id-1)+1];}}}
 if(!threadIdx.x)h.ns[c]=d.ns[c];
}
}
extern "C" int genesis_state_cycle_begin(const double *x,const double *v,double *hx,double *hv,
 const int *counts,int width,int cells,int eligible){
 if(!state_mode()||!eligible){require(!cycle.pending,"host sync during pending maintenance");cycle.retired=false;return 0;}
 require(!cycle.pending,"nested maintenance cycle");
 if(state_mode()!=2&&!production())return 0;
 cycle={true,state_mode()!=2,x,v,width,cells};
 if(state_mode()==3){
  // Migration can activate previously unused slots; check only active outputs later.
  for(size_t i=0;i<3ull*width*cells;++i)hx[i]=hv[i]=std::numeric_limits<double>::quiet_NaN();
  std::puts("GENESIS_STATE_CYCLE: CPU state mirrors poisoned");
 }
 return cycle.retired;
}
extern "C" int genesis_state_cycle_retired(){return cycle.retired;}
extern "C" void genesis_state_cycle_host_refreshed(){cycle.retired=false;}
extern "C" int genesis_state_cycle_consume(const double *x,const double *v,const double *hx,const double *hv,
 const int *counts,int width,int cells,int eligible){
 if(!cycle.pending)return 0;
 require(eligible&&x==cycle.x&&v==cycle.v&&width==cycle.width&&cells==cycle.cells,"invalid retained state handoff");
 ParticleView p;const double *im;
 require(genesis_incoming_metadata_view(&p,&im)&&p.x==x&&p.v==v&&p.width==width&&p.cells==cells,"remapped state storage changed");
 if(state_mode()==2||state_mode()==3){compare_vectors(x,hx,counts,width,cells,state_mode()==3);compare_vectors(v,hv,counts,width,cells,state_mode()==3);
  std::puts(state_mode()==2?"GENESIS_STATE_CYCLE: exact remapped state pass":"GENESIS_STATE_CYCLE: resident state finite and CPU mirrors remain poisoned");}
 cycle.pending=false;return 1;
}
extern "C" int genesis_state_cycle_publish(ParticleView d,ParticleView h,int stage){
 if(!cycle.pending||!cycle.retired)return 0;
 publish_particle_metadata<<<d.cells,128>>>(d,h,stage);checked(cudaGetLastError());return 1;
}
extern "C" int genesis_pack_resident_available(){return pack_mode()!=0;}
extern "C" int genesis_pack_resident_boundary(ParticleView d,ParticleView h,const double *inverse,
 double *host_inverse,double *host_translation){
 if(!pack_mode()||pack_mode()==2||!production())return 0;
 publish_boundary_metadata<<<d.cells,128>>>(d,h,inverse,host_inverse);checked(cudaGetLastError());
 if(pack_mode()==3){for(int c=0;c<d.cells;++c)for(int a=0;a<h.na[c];++a)for(int k=0;k<3;++k)
   host_translation[3ull*(a+d.width*c)+k]=std::numeric_limits<double>::quiet_NaN();
  std::puts("GENESIS_PACK_RESIDENT: CPU translation mirror poisoned");}
 genesis_translation_rebuilt_device();return 1;
}
extern "C" int genesis_pack_resident_copy(const double *x,double *t,double *q,int *n,int *s,
 const double *ht,const double *hq,const int *hn,const int *hs,int width,int cells){
 if(!pack_mode())return 0;
 ParticleView p;const double *translation;const int *counts,*starts;
 if(!genesis_incoming_pair_view(&p,&translation)||!genesis_exclusion_prefix(&counts,&starts))return 0;
 require(p.x==x&&p.width==width&&p.cells==cells,"coordinate pack state layout changed");
 checked(cudaMemcpyAsync(t,translation,24ull*width*cells,cudaMemcpyDeviceToDevice,0));
 checked(cudaMemcpyAsync(q,p.q,8ull*width*cells,cudaMemcpyDeviceToDevice,0));
 checked(cudaMemcpyAsync(n,counts,4ull*cells,cudaMemcpyDeviceToDevice,0));
 checked(cudaMemcpyAsync(s,starts,4ull*cells,cudaMemcpyDeviceToDevice,0));
 if(pack_mode()==2){
  compare_vectors(t,ht,hn,width,cells,false);
  std::vector<double> charge(size_t(width)*cells);std::vector<int> nc(cells),st(cells);
  checked(cudaMemcpy(charge.data(),q,charge.size()*8,cudaMemcpyDeviceToHost));
  checked(cudaMemcpy(nc.data(),n,cells*4,cudaMemcpyDeviceToHost));checked(cudaMemcpy(st.data(),s,cells*4,cudaMemcpyDeviceToHost));
  require(!std::memcmp(nc.data(),hn,cells*4)&&!std::memcmp(st.data(),hs,cells*4),"pack counts/starts differ");
  for(int c=0;c<cells;++c)require(!std::memcmp(charge.data()+size_t(width)*c,hq+size_t(width)*c,8ull*hn[c]),"pack charges differ");
  std::puts("GENESIS_PACK_RESIDENT: exact translation charge counts and starts pass");
 }
 if(pack_mode()==3&&production()){compare_vectors(t,ht,hn,width,cells,true);std::puts("GENESIS_PACK_RESIDENT: GPU translations finite and CPU mirror remains poisoned");}
 return 1;
}
