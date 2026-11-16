#pragma once
#include <cstdint>
#include "precision_types.h"
struct ParticleView {
 genesis::IntegrationReal *x,*v;genesis::WorkingReal *q;genesis::IntegrationReal *m;int *classes,*ids,*solute;std::int16_t *reverse;
 int *ns,*nw,*na,*noh,*groups;genesis::IntegrationReal *distance;int *hlocal,*water;
 int width,cells,water_height,water_capacity,connect,group_height,group_capacity,reverse_capacity,dist_height,dist_capacity;
};
struct AtomRecords {
 const genesis::IntegrationReal *reals;const int *integers,*adds,*exits,*indices;
 int real_rows,real_capacity,int_rows,int_capacity,index_capacity;
};
