#pragma once
#include <cstdint>
struct ParticleView {
 double *x,*v,*q,*m;int *classes,*ids,*solute;std::int16_t *reverse;
 int *ns,*nw,*na,*noh,*groups;double *distance;int *hlocal,*water;
 int width,cells,water_height,water_capacity,connect,group_height,group_capacity,reverse_capacity,dist_height,dist_capacity;
};
struct AtomRecords {
 const double *reals;const int *integers,*adds,*exits,*indices;
 int real_rows,real_capacity,int_rows,int_capacity,index_capacity;
};
