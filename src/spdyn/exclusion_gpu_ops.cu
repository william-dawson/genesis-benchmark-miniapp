#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include "migration_views.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using genesis::WorkingReal;
#define CHECK(c) do{cudaError_t e=(c);if(e!=cudaSuccess){std::fprintf(stderr,"GPU exclusion: %s\n",cudaGetErrorString(e));std::abort();}}while(0)
extern "C" int genesis_bonded_hbm_exclusion(int,const int**,const int**,const signed char**,const int**);
extern "C" int genesis_incoming_pair_view(ParticleView*,const WorkingReal**);
namespace {
struct View {void *p[27];int n[20];};
int4 *mask_records=nullptr;int *mask_counts=nullptr,*inverse=nullptr;int record_cap=0,hash_cap=0;
void *scan_scratch=nullptr;size_t scan_bytes=0;
View gpu{},host{},allocated{};size_t sizes[27]{};bool valid=false;cudaStream_t queue=nullptr;int *error=nullptr,*totals=nullptr;
template<class T> __device__ T* ptr(View v,int i){return static_cast<T*>(v.p[i]);}
__global__ void initialize(View v){
 int c=blockIdx.x*blockDim.x+threadIdx.x;if(c>=v.n[0])return;

 for(int a=0;a<v.n[1];++a)ptr<uint8_t>(v,25)[a+size_t(v.n[1])*(a+v.n[1]*c)]=0;
}
__global__ void atom_total(View v,int *out){if(!threadIdx.x && !blockIdx.x){int c=v.n[0]-1;out[1]=ptr<int>(v,24)[c]+ptr<int>(v,0)[c];}}
struct Count {int ne,nt,nr;__host__ __device__ Count()=default;__host__ __device__ Count(int x):ne(x),nt(x),nr(x){}__host__ __device__ Count(int a,int b,int c):ne(a),nt(b),nr(c){}__host__ __device__ Count operator+(const Count &b)const{return {ne+b.ne,nt+b.nt,nr+b.nr};}};
struct Candidate {int c1,c2,a,b,type,scale,slot,record;};
__device__ Candidate load_global(View v,int c,int first,int last,int type,int scale,double eps,int *bad){
 Candidate q{};q.type=0;if(first<1 || last<1 || first>v.n[19] || last>v.n[19]){atomicExch(bad,5);return q;}
 auto map=ptr<int16_t>(v,3);q.c1=map[2*(first-1)]-1;q.c2=map[2*(last-1)]-1;q.a=map[2*(first-1)+1]-1;q.b=map[2*(last-1)+1]-1;q.type=type;q.scale=scale;
 if(q.c1<0 || q.c1>=v.n[0] || q.c2<0 || q.c2>=v.n[0] || q.a<0 || q.a>=ptr<int>(v,0)[q.c1] || q.b<0 || q.b>=ptr<int>(v,0)[q.c2]){atomicExch(bad,1);q.type=0;return q;}
 if(ptr<int16_t>(v,4)[q.c1+size_t(v.n[5])*q.c2]!=c+1){atomicExch(bad,2);q.type=0;return q;}
 if(type==2 && !(fabs(ptr<WorkingReal>(v,6)[q.a+v.n[7]*q.c1])>eps && fabs(ptr<WorkingReal>(v,6)[q.b+v.n[7]*q.c2])>eps))q.type=4;
 return q;
}
__global__ void construct(View v,double eps,int *bad,int4 *records,int *used,int rcap,int hcap){
 int c=blockIdx.x;if(c>=v.n[0])return;extern __shared__ unsigned long long shared[];
 auto hash=shared;auto ranks=reinterpret_cast<int*>(hash+hcap);
 __shared__ cub::BlockScan<Count,128>::TempStorage scan_storage;
 __shared__ int ne,nt,nr;
 for(int i=threadIdx.x;i<hcap;i+=blockDim.x){hash[i]=0;ranks[i]=2147483647;}
 if(!threadIdx.x){ne=0;nt=0;nr=0;}
 __syncthreads();
 int required=0;
 if(v.n[17]>0){required+=ptr<int>(v,7)[c];if(v.n[18])for(int g=0;g<v.n[11];++g)required+=(g+1)*ptr<int>(v,17)[g+v.n[11]*c];}
 if(v.n[17]>1)required+=3*ptr<int>(v,1)[c]+ptr<int>(v,10)[c];if(v.n[17]>2)required+=ptr<int>(v,13)[c];
 if(required>hcap/2){if(!threadIdx.x)atomicExch(bad,9);return;}
 int phase_offset=0;
 for(int phase=0;phase<5;++phase){
  int count=0;
  if(phase==0 && v.n[17]>0)count=ptr<int>(v,7)[c];
  if(phase==1 && v.n[17]>0 && v.n[18])for(int g=0;g<v.n[11];++g)count+=(g+1)*ptr<int>(v,17)[g+v.n[11]*c];
  if(phase==2 && v.n[17]>1)count=3*ptr<int>(v,1)[c];
  if(phase==3 && v.n[17]>1)count=ptr<int>(v,10)[c];
  if(phase==4 && v.n[17]>2)count=ptr<int>(v,13)[c];
  for(int base=0;base<count;base+=128){int j=base+threadIdx.x;Candidate q{};
   if(j<count){
    if(phase==0 && ptr<int8_t>(v,9)[j+v.n[8]*c]==0){auto p=ptr<int>(v,8)+2*(j+v.n[8]*c);q=load_global(v,c,p[0],p[1],1,0,eps,bad);}
    else if(phase==1){int ordinal=j,g=0;for(;g<v.n[11];++g){int n=(g+1)*ptr<int>(v,17)[g+v.n[11]*c];if(ordinal<n)break;ordinal-=n;}
     if(g>=v.n[11])atomicExch(bad,8);else{int entry=ordinal/(g+1),h=ordinal%(g+1)+1;auto p=ptr<int>(v,18)+v.n[12]*(entry+v.n[13]*(g+v.n[11]*c));q={c,c,p[0]-1,p[h]-1,1,0,0,0};}}
    else if(phase==2){int entry=j/3,pair=j%3;auto p=ptr<int>(v,2)+v.n[3]*(entry+v.n[4]*c);q={c,c,p[pair==2?1:0]-1,p[pair==0?1:2]-1,1,0,0,0};}
    else if(phase==3 && ptr<int8_t>(v,12)[j+v.n[9]*c]==0){auto p=ptr<int>(v,11)+3*(j+v.n[9]*c);q=load_global(v,c,p[0],p[2],2,0,eps,bad);}
    else if(phase==4 && ptr<int8_t>(v,15)[j+v.n[10]*c]==0){auto p=ptr<int>(v,14)+4*(j+v.n[10]*c);q=load_global(v,c,p[0],p[3],3,ptr<int>(v,16)[j+v.n[10]*c]/v.n[16],eps,bad);}
   }
   unsigned bucket=0;int ordinal=phase_offset+j;
   if(q.type){unsigned first=q.c1*v.n[1]+q.a+1,last=q.c2*v.n[1]+q.b+1;
    auto key=(static_cast<unsigned long long>(min(first,last))<<32)|max(first,last);
    bucket=((key>>32)*2654435761u+key)&(hcap-1);int probe=0;
    while(true){auto old=atomicCAS(hash+bucket,0ull,key);if(!old || old==key)break;if(++probe>=hcap){atomicExch(bad,6);q.type=0;break;}bucket=(bucket+1)&(hcap-1);}
    if(q.type)atomicMin(ranks+bucket,ordinal);
   }
   __syncthreads();
   bool accept=q.type && (q.type==1 || (ranks[bucket]==ordinal && !(q.c1==q.c2 && q.a==q.b)));
   Count flags{int(accept && (q.type==1 || q.type==2)),int(accept && q.type==3),int(accept)},offset,total;
   cub::BlockScan<Count,128>(scan_storage).ExclusiveSum(flags,offset,total);
   int exclusion=ne+offset.ne,onefour=nt+offset.nt,record=nr+offset.nr;
   if(accept){int row,col,plane;
    if(q.c1==q.c2){row=max(q.a,q.b);col=min(q.a,q.b);plane=-c-1;}else{plane=ptr<int>(v,5)[q.c1+size_t(v.n[6])*q.c2];row=q.c1<q.c2?q.b:q.a;col=q.c1<q.c2?q.a:q.b;}
    if(plane>v.n[2] || plane==0 || row<0 || col<0 || row>=v.n[1] || col>=v.n[1] || record>=rcap)atomicExch(bad,3);
    else records[size_t(c)*rcap+record]=make_int4(row+1,col+1,plane,0);
   }
   if(flags.ne || flags.nt){int slot=flags.nt?onefour:exclusion,cap=flags.nt?v.n[15]:v.n[14];
    if(slot>=cap)atomicExch(bad,4);else{auto out=ptr<int4>(v,flags.nt?22:21);out[slot+cap*c]=make_int4(q.c1+1,q.c2+1,q.a+1,q.b+1);if(flags.nt)ptr<int>(v,23)[slot+cap*c]=q.scale;}
   }
   __syncthreads();if(!threadIdx.x){ne+=total.ne;nt+=total.nt;nr+=total.nr;}__syncthreads();

  }
  phase_offset+=count;
 }
 if(!threadIdx.x){ptr<int>(v,19)[c]=ne;ptr<int>(v,20)[c]=nt;used[c]=nr;}
}
__global__ void publish(View d,View h){int c=blockIdx.x;
 int ne=min(ptr<int>(d,19)[c],d.n[14]),nt=min(ptr<int>(d,20)[c],d.n[15]);
 if(!threadIdx.x){ptr<int>(h,19)[c]=ne;ptr<int>(h,20)[c]=nt;ptr<int>(h,24)[c]=ptr<int>(d,24)[c];}
 for(int j=threadIdx.x;j<4*ne;j+=blockDim.x)ptr<int>(h,21)[j+4*d.n[14]*c]=ptr<int>(d,21)[j+4*d.n[14]*c];
 for(int j=threadIdx.x;j<4*nt;j+=blockDim.x)ptr<int>(h,22)[j+4*d.n[15]*c]=ptr<int>(d,22)[j+4*d.n[15]*c];
 for(int j=threadIdx.x;j<nt;j+=blockDim.x)ptr<int>(h,23)[j+d.n[15]*c]=ptr<int>(d,23)[j+d.n[15]*c];
}
__global__ void render(View v,const int4 *records,const int *used,int cap){int c=blockIdx.x;
 for(int j=threadIdx.x;j<min(used[c],cap);j+=blockDim.x){int4 r=records[size_t(c)*cap+j];int plane=r.z<0?-r.z-1:r.z-1;
 ptr<uint8_t>(v,r.z<0?25:26)[r.x-1+size_t(v.n[1])*(r.y-1+v.n[1]*plane)]=0;}
}
__global__ void inverse_planes(const int *pairs,int *inv,int cells){int p=blockIdx.x*blockDim.x+threadIdx.x;if(p<cells)inv[pairs[2*p]-1]=p;}
__device__ void clear(uint8_t *packed,uint8_t *dense,size_t i){atomicAnd(reinterpret_cast<unsigned*>(packed)+(i>>5),~(1u<<(i&31)));atomicAnd(reinterpret_cast<unsigned*>(dense)+(i>>2),~(0xffu<<((i&3)*8)));}
__global__ void clear_sparse(uint8_t *packed,uint8_t *dense,const int4 *records,const int *used,int cap,const int *inv,int width,int pitch,int cells){int c=blockIdx.x;
 for(int j=threadIdx.x;j<min(used[c],cap);j+=blockDim.x){int4 r=records[size_t(c)*cap+j];int p=r.z<0?inv[-r.z-1]:cells+r.z-1;
 clear(packed,dense,size_t(p)*pitch+size_t(r.y-1)*width+r.x-1);
 if(r.z<0)clear(packed,dense,size_t(p)*pitch+size_t(r.x-1)*width+r.y-1);}
}
__global__ void diagonals(uint8_t *packed,uint8_t *dense,const int *natom,const int *pairs,int width,int pitch,int cells){int p=blockIdx.x;if(p>=cells)return;
 for(int a=threadIdx.x;a<natom[pairs[2*p]-1];a+=blockDim.x)clear(packed,dense,size_t(p)*pitch+size_t(a)*width+a);
}

}
extern "C" int genesis_exclusion_enabled(){static int mode=[](){const char *s=std::getenv("GENESIS_EXCLUSION_GPU");if(!s || !std::strcmp(s,"0"))return 0;if(!std::strcmp(s,"1"))return 1;if(!std::strcmp(s,"2"))return 2;std::abort();}();return mode;}
extern "C" int genesis_exclusion_active(){return valid;}
extern "C" void genesis_exclusion_build(const View *view,double eps,int *maximum,int *total){
 host=*view;
 if(!queue){CHECK(cudaStreamCreateWithFlags(&queue,cudaStreamNonBlocking));CHECK(cudaMalloc(&error,4));CHECK(cudaMalloc(&totals,8));}
 View fresh=*view;size_t cells=host.n[0],width=host.n[1];
 size_t bytes[27]={4*cells,4*cells,size_t(4)*host.n[3]*host.n[4]*cells,size_t(4)*host.n[19],size_t(2)*host.n[5]*host.n[5],size_t(4)*host.n[6]*host.n[6],sizeof(WorkingReal)*host.n[7]*cells,
 4*cells,size_t(8)*host.n[8]*cells,size_t(host.n[8])*cells,4*cells,size_t(12)*host.n[9]*cells,size_t(host.n[9])*cells,
 4*cells,size_t(16)*host.n[10]*cells,size_t(host.n[10])*cells,size_t(4)*host.n[10]*cells,size_t(4)*host.n[11]*cells,size_t(4)*host.n[12]*host.n[13]*host.n[11]*cells,
 4*cells,4*cells,size_t(16)*host.n[14]*cells,size_t(16)*host.n[15]*cells,size_t(4)*host.n[15]*cells,4*cells,width*width*cells,width*width*host.n[2]};
 void *alias[27]{};ParticleView particle;const WorkingReal *translation;
 if(genesis_incoming_pair_view(&particle,&translation)){
 alias[0]=particle.na;alias[1]=particle.nw;alias[2]=particle.water;alias[3]=particle.reverse;alias[6]=particle.q;alias[17]=particle.hlocal;alias[18]=particle.groups;
 }
 for(int k=1;k<=3;++k){const int *nc,*list,*period;const signed char *types;
 if(genesis_bonded_hbm_exclusion(k,&nc,&list,&types,&period)){int base=k==1?7:(k==2?10:13);alias[base]=const_cast<int*>(nc);alias[base+1]=const_cast<int*>(list);alias[base+2]=const_cast<signed char*>(types);if(k==3)alias[16]=const_cast<int*>(period);}}
 for(int i=0;i<27;++i){if(!allocated.p[i]){sizes[i]=bytes[i];CHECK(cudaMalloc(&allocated.p[i],bytes[i]?bytes[i]:1));}if(sizes[i]!=bytes[i])std::abort();fresh.p[i]=allocated.p[i];
  if(alias[i])fresh.p[i]=alias[i];
  else if(i<19 && (!(i==4 || i==5) || !valid))CHECK(cudaMemcpyAsync(fresh.p[i],host.p[i],bytes[i],cudaMemcpyHostToDevice,queue));
 }
 gpu=fresh;CHECK(cudaMemsetAsync(error,0,4,queue));
 size_t scan=0,reduce=0;auto nc=static_cast<int*>(gpu.p[0]),st=static_cast<int*>(gpu.p[24]);
 CHECK(cub::DeviceScan::ExclusiveSum(nullptr,scan,nc,st,int(cells),queue));CHECK(cub::DeviceReduce::Max(nullptr,reduce,nc,totals,int(cells),queue));
 size_t needed=scan>reduce?scan:reduce;if(needed>scan_bytes){if(scan_scratch)CHECK(cudaFree(scan_scratch));CHECK(cudaMalloc(&scan_scratch,needed));scan_bytes=needed;}
 CHECK(cub::DeviceScan::ExclusiveSum(scan_scratch,scan,nc,st,int(cells),queue));CHECK(cub::DeviceReduce::Max(scan_scratch,reduce,nc,totals,int(cells),queue));
 atom_total<<<1,1,0,queue>>>(gpu,totals);CHECK(cudaGetLastError());
 CHECK(cudaMemsetAsync(gpu.p[19],0,bytes[19],queue));CHECK(cudaMemsetAsync(gpu.p[20],0,bytes[20],queue));
 if(!mask_records){record_cap=host.n[14]+host.n[15]+host.n[9];hash_cap=4096;
 CHECK(cudaMalloc(&mask_records,size_t(cells)*record_cap*sizeof(int4)));CHECK(cudaMalloc(&mask_counts,4*cells));CHECK(cudaMalloc(&inverse,4*cells));
 CHECK(cudaFuncSetAttribute(construct,cudaFuncAttributeMaxDynamicSharedMemorySize,196608));}
 for(;;){CHECK(cudaMemsetAsync(error,0,4,queue));
 construct<<<cells,128,12*hash_cap,queue>>>(gpu,eps,error,mask_records,mask_counts,record_cap,hash_cap);CHECK(cudaGetLastError());
 int code=0;CHECK(cudaMemcpyAsync(&code,error,4,cudaMemcpyDeviceToHost,queue));CHECK(cudaStreamSynchronize(queue));
 if(code==9 && hash_cap<16384){hash_cap*=2;continue;}if(code){std::fprintf(stderr,"GPU exclusion construction error=%d\n",code);std::abort();}break;
 }
 publish<<<cells,128,0,queue>>>(gpu,host);CHECK(cudaGetLastError());
 int bad=0,result[2];CHECK(cudaMemcpyAsync(&bad,error,4,cudaMemcpyDeviceToHost,queue));CHECK(cudaMemcpyAsync(result,totals,8,cudaMemcpyDeviceToHost,queue));

 bool mirror=genesis_exclusion_enabled()==2;const char *verify=std::getenv("GENESIS_EXCLUSION_VERIFY");mirror=mirror || (verify && !std::strcmp(verify,"1"));
 if(mirror){CHECK(cudaMemsetAsync(gpu.p[25],1,bytes[25],queue));CHECK(cudaMemsetAsync(gpu.p[26],1,bytes[26],queue));
 initialize<<<(cells+127)/128,128,0,queue>>>(gpu);CHECK(cudaGetLastError());
 render<<<cells,128,0,queue>>>(gpu,mask_records,mask_counts,record_cap);CHECK(cudaGetLastError());
 for(int i=25;i<27;++i)CHECK(cudaMemcpyAsync(host.p[i],gpu.p[i],bytes[i],cudaMemcpyDeviceToHost,queue));}
 CHECK(cudaStreamSynchronize(queue));if(bad){std::fprintf(stderr,"GPU exclusion validation error=%d\n",bad);std::abort();}
 *maximum=result[0];*total=result[1];valid=true;
}
extern "C" int genesis_exclusion_mask(uint8_t *packed,uint8_t *dense,const int *natom,const int *pairs,int width,int pitch,int local,int near,cudaStream_t stream){
 if(!valid)return 0;if(width>gpu.n[1] || local!=gpu.n[0] || pitch<width*width)std::abort();size_t bytes=size_t(pitch)*near;
 CHECK(cudaMemsetAsync(dense,1,((bytes+7)/8)*8,stream));CHECK(cudaMemsetAsync(packed,0xff,(((bytes+7)/8)+3)&~size_t(3),stream));
 inverse_planes<<<(local+127)/128,128,0,stream>>>(pairs,inverse,local);CHECK(cudaGetLastError());
 diagonals<<<local,128,0,stream>>>(packed,dense,natom,pairs,width,pitch,local);CHECK(cudaGetLastError());
 clear_sparse<<<local,128,0,stream>>>(packed,dense,mask_records,mask_counts,record_cap,inverse,width,pitch,local);CHECK(cudaGetLastError());
 if(bytes%8)CHECK(cudaMemsetAsync(packed+(bytes+7)/8-1,(1u<<(bytes%8))-1,1,stream));return 1;
}
extern "C" int genesis_exclusion_prefix(const int **counts,const int **starts) {
 if(!valid)return 0;*counts=static_cast<const int*>(gpu.p[0]);*starts=static_cast<const int*>(gpu.p[24]);return 1;
}
extern "C" int genesis_exclusion_force_view(int kind,const int **counts,const int **list,const int **starts,int *capacity){
 if(!valid || kind<0 || kind>1)return 0;*counts=static_cast<const int*>(gpu.p[kind?19:20]);*list=static_cast<const int*>(gpu.p[kind?21:22]);*starts=static_cast<const int*>(gpu.p[24]);*capacity=gpu.n[kind?14:15];return 1;
}
