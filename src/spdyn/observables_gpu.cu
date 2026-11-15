#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace {
void require(bool ok,const char *s){if(!ok){std::fprintf(stderr,"GENESIS observables: %s\n",s);std::abort();}}
void checked(cudaError_t e){require(e==cudaSuccess,cudaGetErrorString(e));}
int mode(){static int m=[](){const char *s=std::getenv("GENESIS_OBSERVABLES_GPU");if(!s)return 0;
 require(std::strlen(s)==1&&s[0]>='0'&&s[0]<='3',"mode must be 0, 1, 2 or 3");return s[0]-'0';}();return m;}
bool production(){
 for(const char *name:{"GENESIS_FORCE_VERIFY","GENESIS_STEP_VERIFY","GENESIS_STATE_VERIFY","GENESIS_THERMO_VERIFY","GENESIS_COORD_VERIFY"}){
  const char *s=std::getenv(name);if(s&&std::strcmp(s,"0"))return false;}
 return true;
}
double *parts=nullptr,*total=nullptr,*host=nullptr;int capacity=0,report_step=0;bool pending=false;
__global__ void cell_sums(const double *f,const double *v,const double *m,const int *n,double *out,int width,int cells){
 __shared__ double sums[2][128];int c=blockIdx.x,lane=threadIdx.x;double ff=0,ke=0;
 for(int a=lane;a<n[c];a+=blockDim.x){size_t i=a+size_t(width)*c;double force=0,velocity=0;
  for(int k=0;k<3;++k){force=__dadd_rn(force,__dmul_rn(f[3*i+k],f[3*i+k]));velocity=__dadd_rn(velocity,__dmul_rn(v[3*i+k],v[3*i+k]));}
  ff=__dadd_rn(ff,force);ke=__dadd_rn(ke,__dmul_rn(m[i],velocity));}
 sums[0][lane]=ff;sums[1][lane]=ke;__syncthreads();
 for(int d=64;d;d/=2){if(lane<d)for(int k=0;k<2;++k)sums[k][lane]=__dadd_rn(sums[k][lane],sums[k][lane+d]);__syncthreads();}
 if(!lane){out[2*c]=sums[0][0];out[2*c+1]=sums[1][0];}
}
__global__ void finish_sums(const double *input,double *out,int cells){
 __shared__ double sums[2][128];int lane=threadIdx.x;
 for(int k=0;k<2;++k){double s=0;for(int c=lane;c<cells;c+=blockDim.x)s=__dadd_rn(s,input[2*c+k]);sums[k][lane]=s;}
 __syncthreads();for(int d=64;d;d/=2){if(lane<d)for(int k=0;k<2;++k)sums[k][lane]=__dadd_rn(sums[k][lane],sums[k][lane+d]);__syncthreads();}
 if(!lane){out[0]=sums[0][0];out[1]=sums[1][0];}
}
}
extern "C" int genesis_observables_mode(){return mode();}
extern "C" int genesis_observables_stage(const double *force,const double *vref,const double *mass,
 const int *counts,int width,int cells,int step){
 if(!mode()||(mode()!=2&&!production()))return 0;
 require(!pending&&width>0&&cells>0&&step>0,"invalid reporting stage");
 if(!total){checked(cudaMalloc(&total,16));checked(cudaMallocHost(&host,16));}
 if(capacity<cells){if(parts)checked(cudaFree(parts));checked(cudaMalloc(&parts,16ull*cells));capacity=cells;}
 cell_sums<<<cells,128>>>(force,vref,mass,counts,parts,width,cells);checked(cudaGetLastError());
 finish_sums<<<1,128>>>(parts,total,cells);checked(cudaGetLastError());
 checked(cudaMemcpyAsync(host,total,16,cudaMemcpyDeviceToHost,0));checked(cudaStreamSynchronize(0));
 require(std::isfinite(host[0])&&std::isfinite(host[1])&&host[0]>=0&&host[1]>=0,"invalid report reductions");
 pending=true;report_step=step;return mode()==2?2:1;
}
extern "C" int genesis_observables_take(int step,double *rmsg,double *kinetic){
 if(!pending)return 0;require(step==report_step,"report epoch mismatch");if(mode()==2)return 0;
 *rmsg=host[0];*kinetic=host[1];return 1;
}
extern "C" void genesis_observables_check(int step,double rmsg,double kinetic){
 if(!pending)return;require(step==report_step,"report completion epoch mismatch");
 if(mode()==2){double reference[]={rmsg,kinetic},ratio=0,relative=0;
  for(int k=0;k<2;++k){require(std::isfinite(reference[k]),"nonfinite CPU report");double error=std::abs(host[k]-reference[k]);
   ratio=std::max(ratio,error/(1e-8+1e-10*std::abs(reference[k])));
   relative=std::max(relative,error/std::max(std::abs(reference[k]),1e-300));}
  require(ratio<=1&&relative<=1e-11,"GPU reporting reductions exceed CPU comparison budgets");
  std::printf("GENESIS_OBSERVABLES: CPU reduction reference pass %.12e %.12e\n",ratio,relative);
 }
 pending=false;
}
