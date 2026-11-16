#pragma once
#include <cuda_runtime.h>
#include <cmath>

namespace genesis {
__device__ __forceinline__ float single_periodic_move(float shift,float box) {
    return __fadd_rn(__fmul_rn(box,0.5f),-__fmul_rn(box,roundf(__fdiv_rn(shift,box))));
}
__device__ __forceinline__ int single_cell_coordinate(float x,float origin,float box,float size,int count) {
    float shift=__fadd_rn(x,-origin);
    shift=__fadd_rn(shift,single_periodic_move(shift,box));
    int index=int(__fdiv_rn(shift,size));
    return index==count?index-1:index;
}
__device__ __forceinline__ float single_translation(float x,float origin,float box) {
    return __fadd_rn(single_periodic_move(__fadd_rn(x,-origin),box),-origin);
}
}
