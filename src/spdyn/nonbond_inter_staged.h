#pragma once
#include <cuda_runtime.h>
#include <stdint.h>
#include "precision_types.h"

namespace genesis {
int nonbond_inter_mode();
void launch_nonbond_inter(dim3 grid,dim3 block,size_t shared_bytes,cudaStream_t stream,int mode,

    const WorkingReal *_coord_pbc,          // ( 1:atom_domain, 1:4 )
    WorkingReal *_force,                    // ( 1:atom_domain, 1:3 )
    double *_ene_virial,             // ( 5 )
    double *_ene_viri_mid,           // ( 1:5, 1:univ_maxcell)
    const int8_t *_cell_move,        // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int *_atmcls_pbc,          // ( 1:atom_domain )
    const int *_natom,               // ( 1:ncel_max )
    const int *_start_atom,          // ( 1:ncel_max )
    const WorkingReal *_nonb_lj12,          // ( 1:num_atom_cls, 1:num_atom_cls )
    const WorkingReal *_nonb_lj6,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const WorkingReal *_table_grad,         // ( 1:6*cutoff_int )
    const int *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const int8_t *_univ_mask2,
    const int *_univ_ix_natom,     // ( 1:univ_maxcell1 )
    const uint8_t *_univ_ix_list,  // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int *_univ_iy_natom,     // ( 1:univ_maxcell1 )
    const uint8_t *_univ_iy_list,  // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int *_univ_ij_sort_list, // ( 1:univ_maxcell1 )
    const int8_t *_virial_check,   // ( 1:ncel_max, 1:ncel_max )
    int atom_domain,
    int MaxAtom,
    int MaxAtomCls,
    int num_atom_cls,
    int ncel_local,
    int ncel_bound,
    int ncel_max,
    int cutoff_int,
    int univ_maxcell,
    int univ_maxcell1,
    int univ_ncell_near,
    int univ_mask2_size,
    int univ_natom_max,
    int index_s,
    int index_e,
    int max_iy_natom,
    int check_virial,
    WorkingReal density,
    WorkingReal cutoff2,
    WorkingReal pairlistdist2,
    WorkingReal system_x,
    WorkingReal system_y,
    WorkingReal system_z,
    int inter_tile_log2);
}
