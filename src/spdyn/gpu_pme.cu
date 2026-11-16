#include <cuda_runtime.h>
#include "precision_types.h"
using genesis::WorkingReal;
using genesis::IntegrationReal;
extern "C" int genesis_long_force_requested_();
extern "C" int genesis_long_force_reference();
extern "C" void genesis_long_force_publish(const IntegrationReal *,const int *,const int *,int,int,cudaStream_t);
#include <cstdint>
extern "C" int genesis_coord_copy(void *,int,int,std::uintptr_t,int);
#include <cufft.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>

extern "C" void genesis_pme_mesh_sum_(WorkingReal *, const int *, const int *);
extern "C" void genesis_pme_host_pack_(const int *, const int *, const int *, const int *,
    const int *, const IntegrationReal *, const WorkingReal *, double *, const int *);
extern "C" void genesis_pme_host_scatter_(const int *, const int *, const int *, const int *,
    const int *, const IntegrationReal *, IntegrationReal *, const int *);

namespace {
constexpr int threads = 256;
constexpr double pi = 3.141592653589793238462643383279502884;
#if defined(_SINGLE) || defined(_MIXED)
using MeshComplex = cufftComplex;
using Scale = float3;
constexpr cufftType forward_type = CUFFT_R2C, backward_type = CUFFT_C2R;
#else
using MeshComplex = cufftDoubleComplex;
using Scale = double3;
constexpr cufftType forward_type = CUFFT_D2Z, backward_type = CUFFT_Z2D;
#endif

void check(cudaError_t status, const char *operation) {
  if (status != cudaSuccess) {
    std::fprintf(stderr, "GPU PME: %s: %s\n", operation, cudaGetErrorString(status));
    std::abort();
  }
}
void check(cufftResult status, const char *operation) {
  if (status != CUFFT_SUCCESS) {
    std::fprintf(stderr, "GPU PME: %s: cuFFT status %d\n", operation, int(status));
    std::abort();
  }
}
#define CUDA(call) check((call), #call)
#define FFT(call) check((call), #call)

struct Parameters {
  int nx, ny, nz, ld;
  int x_first, x_count, y_first, y_count, ranks, communicator;
  WorkingReal box[3], alpha, vol4, bs3, bs3d;
#if defined(_SINGLE) || defined(_MIXED)
  WorkingReal scale[3], gfact[3];
#endif
};
struct State {
  Parameters p{};
  bool ready = false;
  bool packed = false;
  bool cached_scale = false;
  bool gather_lanes = false;
  Scale scale{};
  const char *timing_file = nullptr;
  long calls = 0, timed_calls = 0;
  double host_time[8]{}, fortran_time = 0;
  double clear_time = 0;
  int capacity = 0, cells = 0, blocks = 0;
  int host_threads = 0;
  cudaStream_t stream{};
  cufftHandle forward{}, backward{};
  WorkingReal *mesh = nullptr, *theta = nullptr, *b2 = nullptr, *charge = nullptr;
  IntegrationReal *coord = nullptr, *force = nullptr, *host_force = nullptr;
  double *partial = nullptr, *host_partial = nullptr;
  double *packed_input = nullptr, *host_input = nullptr;
  int *natom = nullptr;
  int *host_offsets = nullptr;
  MeshComplex *spectrum = nullptr;
} state;

template<class T> void allocate(T *&pointer, size_t count) {
  CUDA(cudaMalloc(reinterpret_cast<void **>(&pointer), count*sizeof(T)));
}
template<class T> void pinned(T *&pointer, size_t count) {
  CUDA(cudaMallocHost(reinterpret_cast<void **>(&pointer), count*sizeof(T)));
}
template<class T> void release(T *&pointer) {
  if (pointer) CUDA(cudaFree(pointer));
  pointer = nullptr;
}
template<class T> void release_pinned(T *&pointer) {
  if (pointer) CUDA(cudaFreeHost(pointer));
  pointer = nullptr;
}
void release_atoms() {
  release(state.coord); release(state.charge); release(state.force); release(state.natom);
  release(state.packed_input); release_pinned(state.host_input);
  release_pinned(state.host_force);
  release_pinned(state.host_offsets);
  state.capacity = state.cells = 0;
}
void destroy() {
  if (!state.ready) return;
  CUDA(cudaStreamSynchronize(state.stream));
  if (state.timing_file && state.timed_calls) {
    FILE *file=std::fopen(state.timing_file,"w");
    if (!file) { std::perror("PME timing file"); std::abort(); }
    std::fprintf(file,"phase,calls,total_s,mean_us\n");
    const char *labels[]={"pack","h2d_api","kernel_enqueue","d2h_api","synchronize","scatter","scalar_reduce","callback_total"};
    for(int i=0;i<8;++i)
      std::fprintf(file,"%s,%ld,%.17g,%.17g\n",labels[i],state.timed_calls,state.host_time[i],state.host_time[i]*1e6/state.timed_calls);
    std::fprintf(file,"fortran_call,%ld,%.17g,%.17g\n",state.timed_calls,state.fortran_time,state.fortran_time*1e6/state.timed_calls);
    std::fprintf(file,"halo_clear,%ld,%.17g,%.17g\n",state.timed_calls,state.clear_time,state.clear_time*1e6/state.timed_calls);
    if (std::ferror(file) || std::fclose(file)!=0) std::abort();
  }
  FFT(cufftDestroy(state.forward)); FFT(cufftDestroy(state.backward));
  release_atoms();
  release(state.mesh); release(state.theta); release(state.b2); release(state.partial);
  release(state.spectrum); release_pinned(state.host_partial);
  CUDA(cudaStreamDestroy(state.stream));
  state.ready = false;
}

__device__ void wave(int t, Parameters p, int &x, int &y, int &z,
                     WorkingReal &gx, WorkingReal &gy, WorkingReal &gz, WorkingReal &g2) {
  x = t % (p.nx/2+1);
  y = (t/(p.nx/2+1)) % p.ny;
  z = t/((p.nx/2+1)*p.ny);
#if defined(_SINGLE) || defined(_MIXED)
  gx = p.gfact[0]*x;
  gy = p.gfact[1]*(y <= p.ny/2 ? y : y-p.ny);
  gz = p.gfact[2]*(z <= p.nz/2 ? z : z-p.nz);
#else
  gx = 2*pi*x/p.box[0];
  gy = 2*pi*(y <= p.ny/2 ? y : y-p.ny)/p.box[1];
  gz = 2*pi*(z <= p.nz/2 ? z : z-p.nz)/p.box[2];
#endif
  g2 = gx*gx + gy*gy + gz*gz;
}

__global__ void prepare_coefficients(Parameters p, const WorkingReal *b2, WorkingReal *theta, int count) {
  int t = blockIdx.x*blockDim.x + threadIdx.x;
  if (t >= count) return;
  int x,y,z;
  WorkingReal gx,gy,gz,g2;
  wave(t,p,x,y,z,gx,gy,gz,g2);
  WorkingReal exponent = -g2/(4*p.alpha*p.alpha);
  WorkingReal value = 0;
#if defined(_SINGLE) || defined(_MIXED)
  const WorkingReal fact=WorkingReal(0.25)/(-p.alpha*p.alpha);
  exponent=g2*fact;
  const WorkingReal epsilon=1.e-5f;
  if (g2 > epsilon && exponent >= -80.f) {
#else
  if (g2 > 1.e-10 && exponent >= -80.) {
#endif
    value = b2[x]*b2[p.ld+y]*b2[2*p.ld+z]*exp(exponent)/g2;
    if (fabs(value) < 1.e-15) value = 0;
  }
  theta[t] = value;
}

__device__ void spline(WorkingReal u, WorkingReal *w, WorkingReal *d) {
  WorkingReal u2 = u*u, u3 = u2*u, v = 1-u;
  w[0] = u3;
  w[1] = -3*u3 + 3*u2 + 3*u + 1;
  w[2] = 3*u3 - 6*u2 + 4;
  w[3] = v*v*v;
  d[0] = u2;
  d[1] = -3*u2 + 2*u + 1;
  d[2] = 3*u2 - 4*u;
  d[3] = -v*v;
}

template<bool CachedScale, class CoordReal>
__device__ void interpolation(const CoordReal *coord, int atom, Parameters p,
                              int index[3][4], WorkingReal w[3][4], WorkingReal d[3][4], Scale scale) {
  int grid[3] = {p.nx,p.ny,p.nz};
  for (int axis=0; axis<3; ++axis) {
    WorkingReal factor=CachedScale ? (axis==0 ? scale.x : axis==1 ? scale.y : scale.z) : WorkingReal(grid[axis])/p.box[axis];
#if defined(_SINGLE) || defined(_MIXED)
    factor=p.scale[axis];
#endif
    WorkingReal vr = coord[3*atom+axis] * factor;
    // Fortran ANINT rounds halfway cases away from zero.
    vr = vr + grid[axis]/2 - grid[axis]*round(vr/grid[axis]);
    vr = fmin(vr, WorkingReal(grid[axis])-WorkingReal(0.000001));
    int base = min(int(vr),grid[axis]-1);
    spline(vr-base,w[axis],d[axis]);
    for (int k=0; k<4; ++k) {
      int j = base-k;
      if (j < 0) j += grid[axis];
      index[axis][k] = j;
    }
  }
}

template<bool Packed, bool CachedScale, class CoordReal, class ChargeReal>
__global__ void spread(Parameters p, const int *natom, const CoordReal *coord,
                       const ChargeReal *charge, WorkingReal *mesh, int stride, int capacity, Scale scale) {
  int atom = blockIdx.x*blockDim.x + threadIdx.x;
  if (atom >= capacity) return;
  if (!Packed && atom%stride >= natom[atom/stride]) return;
  int index[3][4];
  WorkingReal w[3][4],d[3][4];
  interpolation<CachedScale>(coord,atom,p,index,w,d,scale);
  WorkingReal q = charge[atom];
  for (int z=0; z<4; ++z)
    for (int y=0; y<4; ++y)
      for (int x=0; x<4; ++x) {
        int t = index[0][x] + p.nx*(index[1][y]+p.ny*index[2][z]);
        atomicAdd(mesh+t,w[0][x]*w[1][y]*w[2][z]*q*p.bs3);
      }
}

__global__ void spectral(Parameters p, MeshComplex *spectrum,
                         const WorkingReal *theta, double *partial, int count) {
  __shared__ double sums[4][threads];
  int t = blockIdx.x*blockDim.x + threadIdx.x;
  double e=0,vx=0,vy=0,vz=0;
  if (t < count) {
    int x,y,z;
    WorkingReal gx,gy,gz,g2;
    wave(t,p,x,y,z,gx,gy,gz,g2);
    auto value = spectrum[t];
#if defined(_SINGLE) || defined(_MIXED)
    WorkingReal fact=WorkingReal(0.25)/(-p.alpha*p.alpha);
    WorkingReal factor=g2 > 1.e-5f ? -WorkingReal(2)*(WorkingReal(1)-g2*fact)/g2 : 0;
#else
    double factor = g2 > 1.e-10 ? -2*(1+g2/(4*p.alpha*p.alpha))/g2 : 0;
#endif
    double weight = (x==0 || x==p.nx/2) ? p.vol4*0.5 : p.vol4;
    bool owned_x = (x>=p.x_first && x<p.x_first+p.x_count) || (x==p.nx/2 && p.x_first==0);
    bool owned_y = y>=p.y_first && y<p.y_first+p.y_count;
    if (!owned_x || !owned_y) weight=0;
    e = (value.x*value.x + value.y*value.y)*theta[t]*weight;
    vx=e*(1+factor*gx*gx); vy=e*(1+factor*gy*gy); vz=e*(1+factor*gz*gz);
    // CPU CMPLX(theta) omits KIND and narrows the force filter to float.
    // Preserve that operation; energy and virial still use double theta.
    WorkingReal force_theta = WorkingReal(float(theta[t]));
#if defined(_SINGLE) || defined(_MIXED)
    spectrum[t] = make_cuFloatComplex(value.x*force_theta,value.y*force_theta);
#else
    spectrum[t] = make_cuDoubleComplex(value.x*force_theta,value.y*force_theta);
#endif
  }
  sums[0][threadIdx.x]=e; sums[1][threadIdx.x]=vx;
  sums[2][threadIdx.x]=vy; sums[3][threadIdx.x]=vz;
  __syncthreads();
  for (int offset=threads/2; offset; offset/=2) {
    if (threadIdx.x < offset)
      for (int k=0;k<4;++k) sums[k][threadIdx.x]+=sums[k][threadIdx.x+offset];
    __syncthreads();
  }
  if (threadIdx.x==0)
    for (int k=0;k<4;++k) partial[4*blockIdx.x+k]=sums[k][0];
}

template<bool Packed, bool CachedScale, class CoordReal, class ChargeReal>
__global__ void gather(Parameters p, const int *natom, const CoordReal *coord,
                       const ChargeReal *charge, const WorkingReal *mesh, IntegrationReal *force,
                       int stride, int capacity, Scale scale) {
  int atom = blockIdx.x*blockDim.x + threadIdx.x;
  if (atom >= capacity) return;
  if (!Packed && atom%stride >= natom[atom/stride]) return;
  int index[3][4];
  WorkingReal w[3][4],d[3][4], f[3] = {0,0,0};
  interpolation<CachedScale>(coord,atom,p,index,w,d,scale);
  for (int z=0;z<4;++z)
    for (int y=0;y<4;++y)
      for (int x=0;x<4;++x) {
        WorkingReal v=mesh[index[0][x]+p.nx*(index[1][y]+p.ny*index[2][z])];
        f[0]+=d[0][x]*w[1][y]*w[2][z]*v;
        f[1]+=w[0][x]*d[1][y]*w[2][z]*v;
        f[2]+=w[0][x]*w[1][y]*d[2][z]*v;
      }
  int grid[3] = {p.nx,p.ny,p.nz};
  for (int k=0;k<3;++k) {
    WorkingReal factor=CachedScale ? (k==0 ? scale.x : k==1 ? scale.y : scale.z) : WorkingReal(grid[k])/p.box[k];
#if defined(_SINGLE) || defined(_MIXED)
    factor=p.scale[k];
#endif
    force[3*atom+k] = -f[k]*WorkingReal(charge[atom])*factor*p.vol4*p.bs3d;
  }
}

template<bool Packed, class CoordReal, class ChargeReal> void launch_spread(State &s,const CoordReal *coord,const ChargeReal *charge,int stride,int atoms) {
  if (s.cached_scale)
    spread<Packed,true><<<(atoms+threads-1)/threads,threads,0,s.stream>>>(s.p,s.natom,coord,charge,s.mesh,stride,atoms,s.scale);
  else
    spread<Packed,false><<<(atoms+threads-1)/threads,threads,0,s.stream>>>(s.p,s.natom,coord,charge,s.mesh,stride,atoms,s.scale);
}

template<class CoordReal, class ChargeReal>
__global__ void gather_four(Parameters p, const CoordReal *coord, const ChargeReal *charge,
                            const WorkingReal *mesh, IntegrationReal *force, int atoms, Scale scale) {
  int thread=blockIdx.x*blockDim.x+threadIdx.x;
  int atom=thread/4, lane=thread%4;
  if (atom>=atoms) return;
  int index[3][4];
  WorkingReal w[3][4],d[3][4];
  interpolation<true>(coord,atom,p,index,w,d,scale);
  // Four lanes share one atom; adjacent lanes read adjacent mesh x points.
  int x=lane==0 ? index[0][0] : lane==1 ? index[0][1] : lane==2 ? index[0][2] : index[0][3];
  WorkingReal wx=lane==0 ? w[0][0] : lane==1 ? w[0][1] : lane==2 ? w[0][2] : w[0][3];
  WorkingReal dx=lane==0 ? d[0][0] : lane==1 ? d[0][1] : lane==2 ? d[0][2] : d[0][3];
  WorkingReal fx=0,fy=0,fz=0;
  for(int z=0;z<4;++z)
    for(int y=0;y<4;++y) {
      WorkingReal v=mesh[x+p.nx*(index[1][y]+p.ny*index[2][z])];
      fx+=dx*w[1][y]*w[2][z]*v;
      fy+=wx*d[1][y]*w[2][z]*v;
      fz+=wx*w[1][y]*d[2][z]*v;
    }
  unsigned mask=__activemask();
  for(int offset=2;offset;offset/=2) {
    fx+=__shfl_down_sync(mask,fx,offset,4);
    fy+=__shfl_down_sync(mask,fy,offset,4);
    fz+=__shfl_down_sync(mask,fz,offset,4);
  }
  if(lane==0) {
    force[3*atom]=-fx*WorkingReal(charge[atom])*scale.x*p.vol4*p.bs3d;
    force[3*atom+1]=-fy*WorkingReal(charge[atom])*scale.y*p.vol4*p.bs3d;
    force[3*atom+2]=-fz*WorkingReal(charge[atom])*scale.z*p.vol4*p.bs3d;
  }
}

template<bool Packed, class CoordReal, class ChargeReal> void launch_gather(State &s,const CoordReal *coord,const ChargeReal *charge,int stride,int atoms) {
  if (s.gather_lanes)
    gather_four<<<(4*atoms+threads-1)/threads,threads,0,s.stream>>>(s.p,coord,charge,s.mesh,s.force,atoms,s.scale);
  else if (s.cached_scale)
    gather<Packed,true><<<(atoms+threads-1)/threads,threads,0,s.stream>>>(s.p,s.natom,coord,charge,s.mesh,s.force,stride,atoms,s.scale);
  else
    gather<Packed,false><<<(atoms+threads-1)/threads,threads,0,s.stream>>>(s.p,s.natom,coord,charge,s.mesh,s.force,stride,atoms,s.scale);
}
double wall_time() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
} // namespace

extern "C" int genesis_pme_gpu_raw_() {
  const char *value=std::getenv("GENESIS_PME_RAW");
  return value && std::strcmp(value,"1")==0;
}
extern "C" double genesis_pme_gpu_wall_time_() { return state.timing_file ? wall_time() : 0.; }
extern "C" int genesis_pme_gpu_parallel_clear_() {
  const char *value=std::getenv("GENESIS_PME_CLEAR");
  return value && std::strcmp(value,"1")==0;
}
extern "C" void genesis_pme_gpu_clear_time_(const double *elapsed) {
  if (state.timing_file && state.calls>=60) state.clear_time+=*elapsed;
}
extern "C" void genesis_pme_gpu_fortran_time_(const double *elapsed) {
  if (state.timing_file && state.calls>60) state.fortran_time+=*elapsed;
}

extern "C" int genesis_pme_gpu_enabled_() {
  static int enabled = [] {
    const char *value=std::getenv("GENESIS_PME_GPU");
    return value && std::strcmp(value,"1")==0 ? 1 : 0;
  }();
  return enabled;
}

extern "C" void genesis_pme_gpu_prepare_(const int *grid, const int *ld, const WorkingReal *b2,
    const WorkingReal *box, const WorkingReal *alpha, const WorkingReal *el_fact,
    const WorkingReal *bs3, const WorkingReal *bs3d, const int *xy_rank, const int *xy_size,
    const int *z_rank, const int *z_size, const int *ranks, const int *communicator,
    const WorkingReal *prepared_scale, const WorkingReal *box_inv, const WorkingReal *vol4) {
  auto &s=state;
  if (s.ready && (s.p.nx!=grid[0] || s.p.ny!=grid[1] || s.p.nz!=grid[2] || s.p.ld!=*ld)) destroy();
  if (!s.ready) {
    s.p.nx=grid[0]; s.p.ny=grid[1]; s.p.nz=grid[2]; s.p.ld=*ld;
    size_t real_count=size_t(grid[0])*grid[1]*grid[2];
    size_t complex_count=size_t(grid[0]/2+1)*grid[1]*grid[2];
    s.blocks=(complex_count+threads-1)/threads;
    const char *packed=std::getenv("GENESIS_PME_PACKED");
    s.packed=packed && std::strcmp(packed,"1")==0;
    const char *scale=std::getenv("GENESIS_PME_SCALE");
    s.cached_scale=scale && std::strcmp(scale,"1")==0;
    const char *lanes=std::getenv("GENESIS_PME_GATHER_LANES");
    s.gather_lanes=lanes && std::strcmp(lanes,"4")==0;
    if (s.gather_lanes && (!s.packed || !s.cached_scale)) {
      std::fprintf(stderr,"GPU PME: four-lane gather requires packed atoms and cached scales\n");
      std::abort();
    }
    s.timing_file=std::getenv("GENESIS_PME_TIMING_FILE");
    const char *workers=std::getenv("GENESIS_PME_HOST_THREADS");
    if (workers) {
      char *end=nullptr;
      long value=std::strtol(workers,&end,10);
      if (!*workers || *end || value<0 || value>32) {
        std::fprintf(stderr,"GPU PME: host threads must be between 0 and 32\n");
        std::abort();
      }
      s.host_threads=int(value);
    }
    int priority=0, lowest=0, highest=0;
    const char *urgent=std::getenv("GENESIS_PME_PRIORITY");
    if (urgent && std::strcmp(urgent,"1")==0) {
      CUDA(cudaDeviceGetStreamPriorityRange(&lowest,&highest));
      priority=highest;
    }
    CUDA(cudaStreamCreateWithPriority(&s.stream,cudaStreamNonBlocking,priority));
    allocate(s.mesh,real_count); allocate(s.spectrum,complex_count);
    allocate(s.theta,complex_count); allocate(s.b2,3*size_t(*ld));
    allocate(s.partial,4*size_t(s.blocks)); pinned(s.host_partial,4*size_t(s.blocks));
    FFT(cufftPlan3d(&s.forward,grid[2],grid[1],grid[0],forward_type));
    FFT(cufftPlan3d(&s.backward,grid[2],grid[1],grid[0],backward_type));
    FFT(cufftSetStream(s.forward,s.stream)); FFT(cufftSetStream(s.backward,s.stream));
    CUDA(cudaMemcpyAsync(s.b2,b2,3*size_t(*ld)*sizeof(WorkingReal),cudaMemcpyHostToDevice,s.stream));
    s.ready=true;
    std::fprintf(stderr,"GENESIS_GPU_PME: resident order-4 %zu-bit mesh %d %d %d; reusable cuFFT plans\n",8*sizeof(WorkingReal),grid[0],grid[1],grid[2]);
  }
  for(int k=0;k<3;++k) s.p.box[k]=box[k];
#if defined(_SINGLE) || defined(_MIXED)
  s.scale=make_float3(prepared_scale[0],prepared_scale[1],prepared_scale[2]);
  for(int k=0;k<3;++k) {
    s.p.scale[k]=prepared_scale[k];
    s.p.gfact[k]=WorkingReal(2)*WorkingReal(pi)*box_inv[k];
  }
#else
  s.scale=make_double3(double(grid[0])/box[0],double(grid[1])/box[1],double(grid[2])/box[2]);
#endif
  s.p.alpha=*alpha; s.p.bs3=*bs3; s.p.bs3d=*bs3d;
  s.p.x_count=(grid[0]/2)/(*xy_size); s.p.x_first=(*xy_rank)*s.p.x_count;
  s.p.y_count=grid[1]/(*z_size); s.p.y_first=(*z_rank)*s.p.y_count;
  s.p.ranks=*ranks; s.p.communicator=*communicator;
#if defined(_SINGLE) || defined(_MIXED)
  s.p.vol4=*vol4;
#else
  s.p.vol4=4*pi*(*el_fact)/(box[0]*box[1]*box[2]);
#endif
  int count=(grid[0]/2+1)*grid[1]*grid[2];
  prepare_coefficients<<<s.blocks,threads,0,s.stream>>>(s.p,s.b2,s.theta,count);
  CUDA(cudaGetLastError());
}

extern "C" void genesis_pme_gpu_recip_(const int *ncell, const int *stride,
    const int *natom, const IntegrationReal *coord, const WorkingReal *charge, IntegrationReal *force, double *result) {
  auto &s=state;
  bool measure=s.timing_file && s.calls>=60;
  double times[8]{};
  if (measure) times[0]=wall_time();
  if (!s.ready) { std::fprintf(stderr,"GPU PME: prepare was not called\n"); std::abort(); }
  int capacity=(*ncell)*(*stride);
  if (capacity!=s.capacity || *ncell!=s.cells) {
    CUDA(cudaStreamSynchronize(s.stream));
    release_atoms();
    s.capacity=capacity; s.cells=*ncell;
    if (s.packed) {
      allocate(s.packed_input,4*size_t(capacity)); pinned(s.host_input,4*size_t(capacity));
      pinned(s.host_offsets,size_t(*ncell)+1);
    } else {
      allocate(s.coord,3*size_t(capacity)); allocate(s.charge,size_t(capacity));
      allocate(s.natom,size_t(*ncell));
    }
    allocate(s.force,3*size_t(capacity));
    pinned(s.host_force,3*size_t(capacity));
  }
  int atoms=capacity;
  bool resident_input=false;
  if (s.packed) {
    s.host_offsets[0]=0;
    for(int cell=0;cell<*ncell;++cell) s.host_offsets[cell+1]=s.host_offsets[cell]+natom[cell];
    atoms=s.host_offsets[*ncell];
    resident_input=genesis_coord_copy(s.packed_input,atoms,4,reinterpret_cast<std::uintptr_t>(s.stream),2)!=0;
    if(!resident_input) {
    if (s.host_threads) {
      genesis_pme_host_pack_(ncell,stride,&atoms,natom,s.host_offsets,coord,charge,s.host_input,&s.host_threads);
    } else {
    int offset=0;
    for(int cell=0;cell<*ncell;++cell) {
#if defined(_SINGLE) || defined(_MIXED)
      for(int atom=0;atom<natom[cell];++atom) {
        for(int k=0;k<3;++k) s.host_input[3*(offset+atom)+k]=coord[3*(cell*(*stride)+atom)+k];
        s.host_input[3*atoms+offset+atom]=charge[cell*(*stride)+atom];
      }
#else
      std::memcpy(s.host_input+3*offset,coord+3*cell*(*stride),3*size_t(natom[cell])*sizeof(double));
      std::memcpy(s.host_input+3*atoms+offset,charge+cell*(*stride),size_t(natom[cell])*sizeof(double));
#endif
      offset+=natom[cell];
    }
    }
    }
  }
  if (measure) times[1]=wall_time();
  if (s.packed) {
    if(!resident_input) CUDA(cudaMemcpyAsync(s.packed_input,s.host_input,4*size_t(atoms)*sizeof(double),cudaMemcpyHostToDevice,s.stream));
  } else {
    CUDA(cudaMemcpyAsync(s.coord,coord,3*size_t(capacity)*sizeof(IntegrationReal),cudaMemcpyHostToDevice,s.stream));
    CUDA(cudaMemcpyAsync(s.charge,charge,size_t(capacity)*sizeof(WorkingReal),cudaMemcpyHostToDevice,s.stream));
    CUDA(cudaMemcpyAsync(s.natom,natom,size_t(*ncell)*sizeof(int),cudaMemcpyHostToDevice,s.stream));
  }
  if (measure) times[2]=wall_time();
  CUDA(cudaMemsetAsync(s.mesh,0,size_t(s.p.nx)*s.p.ny*s.p.nz*sizeof(WorkingReal),s.stream));
  if (s.packed)
    launch_spread<true>(s,s.packed_input,s.packed_input+3*atoms,*stride,atoms);
  else
    launch_spread<false>(s,s.coord,s.charge,*stride,atoms);
  CUDA(cudaGetLastError());
  if (s.p.ranks>1) {
    // Each rank spreads owned charges once; keep the reduced mesh on-device.
    CUDA(cudaStreamSynchronize(s.stream));
    int mesh_count=s.p.nx*s.p.ny*s.p.nz;
    genesis_pme_mesh_sum_(s.mesh,&mesh_count,&s.p.communicator);
  }
#if defined(_SINGLE) || defined(_MIXED)
  FFT(cufftExecR2C(s.forward,s.mesh,s.spectrum));
#else
  FFT(cufftExecD2Z(s.forward,s.mesh,s.spectrum));
#endif
  int count=(s.p.nx/2+1)*s.p.ny*s.p.nz;
  spectral<<<s.blocks,threads,0,s.stream>>>(s.p,s.spectrum,s.theta,s.partial,count);
  CUDA(cudaGetLastError());
#if defined(_SINGLE) || defined(_MIXED)
  FFT(cufftExecC2R(s.backward,s.spectrum,s.mesh));
#else
  FFT(cufftExecZ2D(s.backward,s.spectrum,s.mesh));
#endif
  if (s.packed)
    launch_gather<true>(s,s.packed_input,s.packed_input+3*atoms,*stride,atoms);
  else
    launch_gather<false>(s,s.coord,s.charge,*stride,atoms);
  CUDA(cudaGetLastError());
  if (measure) times[3]=wall_time();
  const bool resident_force=genesis_long_force_requested_()!=0;
  if(resident_force) {
    if(!s.packed||s.p.ranks!=1) {
      std::fprintf(stderr,"Resident long force requires one-rank packed PME\n");std::abort();
    }
    genesis_long_force_publish(s.force,natom,s.host_offsets,*stride,*ncell,s.stream);
  }
  const bool host_force=!resident_force||genesis_long_force_reference();
  if(host_force) CUDA(cudaMemcpyAsync(s.host_force,s.force,3*size_t(atoms)*sizeof(IntegrationReal),cudaMemcpyDeviceToHost,s.stream));
  CUDA(cudaMemcpyAsync(s.host_partial,s.partial,4*size_t(s.blocks)*sizeof(double),cudaMemcpyDeviceToHost,s.stream));
  if (measure) times[4]=wall_time();
  CUDA(cudaStreamSynchronize(s.stream));
  if (measure) times[5]=wall_time();
  if(host_force) {
  if (s.packed && s.host_threads) {
    genesis_pme_host_scatter_(ncell,stride,&atoms,natom,s.host_offsets,s.host_force,force,&s.host_threads);
  } else {
  int offset=0;
  for(int cell=0;cell<*ncell;++cell) {
    std::memcpy(force+3*cell*(*stride),s.host_force+3*(s.packed ? offset : cell*(*stride)),3*size_t(natom[cell])*sizeof(IntegrationReal));
    offset+=natom[cell];
  }
  }
  }
  if (measure) times[6]=wall_time();
  for(int k=0;k<4;++k) result[k]=0;
  for(int block=0;block<s.blocks;++block)
    for(int k=0;k<4;++k) result[k]+=s.host_partial[4*block+k];
  if (measure) {
    times[7]=wall_time();
    for(int i=0;i<7;++i) s.host_time[i]+=times[i+1]-times[i];
    s.host_time[7]+=times[7]-times[0];
    ++s.timed_calls;
  }
  ++s.calls;
}

extern "C" void genesis_pme_gpu_release_() { destroy(); }
