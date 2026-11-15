#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "migration_views.h"

extern "C" int genesis_incoming_metadata_view(ParticleView*,const double**);
namespace {
void require(bool ok,const char *message) {
    if(!ok){std::fprintf(stderr,"GENESIS metadata import: %s\n",message);std::abort();}
}
void checked(cudaError_t code){require(code==cudaSuccess,cudaGetErrorString(code));}
int mode(){static int value=[] {
    const char *s=std::getenv("GENESIS_METADATA_IMPORT");if(!s)return 0;
    require(std::strlen(s)==1&&s[0]>='0'&&s[0]<='3',"mode must be 0, 1, 2 or 3");return s[0]-'0';
}();return value;}
void compare(const std::vector<unsigned char> &actual,const void *reference,size_t offset,size_t bytes) {
    require(std::memcmp(actual.data()+offset,static_cast<const unsigned char*>(reference)+offset,bytes)==0,
            "active field differs from CPU metadata");
}
}
extern "C" int genesis_metadata_import(double *mass,double *inverse,int *counts,int *waters,int *water_list,
    int *solute,int *group_counts,int *groups,double *distances,
    const double *hmass,const double *hinverse,const int *hcounts,const int *hwaters,const int *hwater_list,
    const int *hsolute,const int *hgroup_counts,const int *hgroups,const double *hdistances,
    int width,int cells,int wh,int wc,int connect,int gh,int gc,int dh,int dc,int eligible) {
    if(!mode()||!eligible)return 0;
    ParticleView p;const double *im=nullptr;
    require(genesis_incoming_metadata_view(&p,&im),"resident metadata missing after maintenance");
    require(p.width==width&&p.cells==cells&&p.water_height==wh&&p.water_capacity==wc&&p.connect==connect&&
            p.group_height==gh&&p.group_capacity==gc&&p.dist_height==dh&&p.dist_capacity==dc,
            "resident metadata shape differs from integration");
    void *out[]={mass,inverse,counts,waters,water_list,solute,group_counts,groups,distances};
    const void *in[]={p.m,im,p.na,p.nw,p.water,p.noh,p.hlocal,p.groups,p.distance};
    const void *host[]={hmass,hinverse,hcounts,hwaters,hwater_list,hsolute,hgroup_counts,hgroups,hdistances};
    size_t sizes[]={8ull*width*cells,8ull*width*cells,4ull*cells,4ull*cells,4ull*wh*wc*cells,
                    4ull*cells,4ull*connect*cells,4ull*gh*gc*connect*cells,8ull*dh*dc*connect*cells};
    for(int field=0;field<9;++field) {
        if(mode()==3)checked(cudaMemsetAsync(out[field],0xff,sizes[field],0));
        checked(cudaMemcpyAsync(out[field],in[field],sizes[field],cudaMemcpyDeviceToDevice,0));
    }
    if(mode()==2) {
        for(int field=0;field<9;++field) {
            std::vector<unsigned char> actual(sizes[field]);
            checked(cudaMemcpyAsync(actual.data(),out[field],sizes[field],cudaMemcpyDeviceToHost,0));
            checked(cudaStreamSynchronize(0));
            if(field==2||field==3||field==5||field==6)compare(actual,host[field],0,sizes[field]);
            else for(int c=0;c<cells;++c) {
                require(hcounts[c]>=0&&hcounts[c]<=width&&hwaters[c]>=0&&hwaters[c]<=wc,"invalid CPU count");
                if(field<2)compare(actual,host[field],8ull*width*c,8ull*hcounts[c]);
                else if(field==4)compare(actual,host[field],4ull*wh*wc*c,4ull*wh*hwaters[c]);
                else for(int g=0;g<connect;++g) {
                    int n=hgroup_counts[g+connect*c];
                    require(n>=0&&n<=gc&&n<=dc,"invalid CPU group count");
                    const int height=field==7?gh:dh,capacity=field==7?gc:dc,bytes=field==7?4:8;
                    // Padding sites need not be initialized; compare only each group's live sites.
                    for(int j=0;j<n;++j)compare(actual,host[field],size_t(bytes)*height*(j+capacity*(g+connect*c)),
                                               size_t(bytes)*(g+2));
                }
            }
        }
        std::puts("GENESIS_METADATA_IMPORT: exact active metadata pass");
    }
    if(mode()==3)std::puts("GENESIS_METADATA_IMPORT: poisoned device metadata overwritten");
    return 1;
}
