#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <limits>
#include "precision_types.h"
using genesis::IntegrationReal;

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
#if defined(_SINGLE)
double report_relative_bound[2]{};
#endif
__global__ void cell_sums(const IntegrationReal *f,const IntegrationReal *v,const IntegrationReal *m,const int *n,double *out,int width,int cells){
 __shared__ double sums[2][128];int c=blockIdx.x,lane=threadIdx.x;double ff=0,ke=0;
 for(int a=lane;a<n[c];a+=blockDim.x){size_t i=a+size_t(width)*c;
#if defined(_SINGLE)
  double force=0;float velocity=0;
  for(int k=0;k<3;++k){force=__dadd_rn(force,__fmul_rn(f[3*i+k],f[3*i+k]));velocity=__fadd_rn(velocity,__fmul_rn(v[3*i+k],v[3*i+k]));}
  ff=__dadd_rn(ff,force);ke=__dadd_rn(ke,__fmul_rn(m[i],velocity));}
#else
  double force=0,velocity=0;
  for(int k=0;k<3;++k){force=__dadd_rn(force,__dmul_rn(f[3*i+k],f[3*i+k]));velocity=__dadd_rn(velocity,__dmul_rn(v[3*i+k],v[3*i+k]));}
  ff=__dadd_rn(ff,force);ke=__dadd_rn(ke,__dmul_rn(m[i],velocity));}
#endif
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
extern "C" int genesis_observables_stage(const IntegrationReal *force,const IntegrationReal *vref,const IntegrationReal *mass,
 const int *counts,int width,int cells,int step){
 if(!mode()||(mode()!=2&&!production()))return 0;
 require(!pending&&width>0&&cells>0&&step>0,"invalid reporting stage");
#if defined(_SINGLE)
 const double u=0.5*std::numeric_limits<float>::epsilon();
 const double ud=0.5*std::numeric_limits<double>::epsilon();
 const double work=6.0*width*cells+2048.0;
 require(work*ud<1,"report reduction too large for rounding bound");
 const double gd=work*ud/(1-work*ud);
 for(int k=0;k<2;++k){const double operations=k?5.0:1.0;
  const double gf=operations*u/(1-operations*u);
  const double beta=gf+gd+gf*gd;
  report_relative_bound[k]=2*beta/(1-beta);}
#endif
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
#if defined(_SINGLE)
   ratio=std::max(ratio,error/(1e-8+report_relative_bound[k]*std::abs(reference[k])));
   require(error<=report_relative_bound[k]*std::abs(reference[k]),"GPU report exceeds FP32 rounding envelope");
#else
   ratio=std::max(ratio,error/(1e-8+1e-10*std::abs(reference[k])));
#endif
   relative=std::max(relative,error/std::max(std::abs(reference[k]),1e-300));}
#if defined(_SINGLE)
  require(ratio<=1,"GPU reporting reductions exceed FP32 comparison budgets");
#else
  require(ratio<=1&&relative<=1e-11,"GPU reporting reductions exceed CPU comparison budgets");
#endif
  std::printf("GENESIS_OBSERVABLES: CPU reduction reference pass %.12e %.12e\n",ratio,relative);
 }
 pending=false;
}
