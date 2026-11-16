#include <cstdio>
#include <cstdlib>
#include "gpu_common.h"
#include "nonbond_inter_staged.h"
using REAL = genesis::WorkingReal;
extern __shared__ REAL smem[];

#define coord_pbc(X,Y)            _coord_pbc           [CALIDX2((X)-1,atom_domain, (Y)-1,4)]
#define force(X,Y)                _force               [CALIDX2((X)-1,atom_domain, (Y)-1,3)]
#define atmcls_pbc(X)             _atmcls_pbc          [(X)-1]
#define cell_move(X,Y,Z)          _cell_move           [CALIDX3((X)-1,3, (Y)-1,ncel_max, (Z)-1,ncel_max)]
#define natom(Z)                  _natom               [(Z)-1]
#define start_atom(Z)             _start_atom          [(Z)-1]
#define nonb_lj12(Y,Z)            _nonb_lj12           [CALIDX2((Y)-1,num_atom_cls, (Z)-1,num_atom_cls)]
#define nonb_lj6(Y,Z)             _nonb_lj6            [CALIDX2((Y)-1,num_atom_cls, (Z)-1,num_atom_cls)]
#define nonb_lj6_factor(Z)        _nonb_lj6_factor     [(Z)-1]
#define table_ene(Z)              _table_ene           [(Z)-1]
#define table_grad(Z)             _table_grad          [(Z)-1]
#define univ_cell_pairlist1(Y,Z)  _univ_cell_pairlist1 [CALIDX2((Y)-1,2, (Z)-1,univ_maxcell)]
#define univ_mask2(Y,Z)           _univ_mask2          [CALIDX2((Y)-1,univ_mask2_size, (Z)-1,univ_ncell_near)]
#define univ_ix_list(Y,Z)         _univ_ix_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define univ_iy_list(Y,Z)         _univ_iy_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define univ_ix_natom(Z)          _univ_ix_natom       [(Z)-1]
#define univ_iy_natom(Z)          _univ_iy_natom       [(Z)-1]
#define univ_ij_sort_list(Z)      _univ_ij_sort_list   [(Z)-1]
#define virial_check(Y,Z)         _virial_check        [CALIDX2((Y)-1,ncel_max, (Z)-1,ncel_max)]

#define ene_virial(Z)      _ene_virial[(Z)-1]
#define ene_viri_mid(Y,Z)  _ene_viri_mid[CALIDX2((Y)-1,5, (Z)-1,univ_maxcell)]
#define force_local(Z)     _force_local[(Z)-1]

__launch_bounds__(128,7)
__global__ void nonbond_inter_staged_force(
    const REAL *_coord_pbc,          // ( 1:atom_domain, 1:4 )
    REAL *_force,                    // ( 1:atom_domain, 1:3 )
    double *_ene_virial,             // ( 5 )
    double *_ene_viri_mid,           // ( 1:5, 1:univ_maxcell)
    const int8_t *_cell_move,        // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int *_atmcls_pbc,          // ( 1:atom_domain )
    const int *_natom,               // ( 1:ncel_max )
    const int *_start_atom,          // ( 1:ncel_max )
    const REAL *_nonb_lj12,          // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL *_nonb_lj6,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL *_table_grad,         // ( 1:6*cutoff_int )
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
    REAL density,
    REAL cutoff2,
    REAL pairlistdist2,
    REAL system_x,
    REAL system_y,
    REAL system_z,
    int inter_tile_log2)
{
#undef num_thread_xx
#undef num_thread_xy
#undef id_thread_xx
#undef id_thread_xy
#define num_thread_xx 8
#define num_thread_xy 4
#define id_thread_xx (id_thread_x / num_thread_xy)
#define id_thread_xy (id_thread_x & (num_thread_xy - 1))

    const int num_thread_x = blockDim.x;
    const int num_thread_y = blockDim.y;
    const int id_thread_x = threadIdx.x;
    const int id_thread_y = threadIdx.y;
    REAL _force_local[3];

    /* shared memory */
    REAL *_force_iy_smem = &smem[id_thread_y * (max_iy_natom * 7)];
    REAL *j_values = _force_iy_smem + 3 * max_iy_natom;
    int *j_ids = reinterpret_cast<int *>(smem + 7 * num_thread_y * max_iy_natom)
                 + id_thread_y * max_iy_natom;
    int *j_classes = j_ids + num_thread_y * max_iy_natom;
#define j_value(D) j_values[((D)-1)*max_iy_natom+iiy-iiy_s]
#define force_iy_smem(X, Y) _force_iy_smem[CALIDX2((X) - 1, 3, (Y) - 1, max_iy_natom)]

    const int tile = blockIdx.x & ((1 << inter_tile_log2) - 1);
    int index = index_s + id_thread_y + num_thread_y * (blockIdx.x >> inter_tile_log2);
    if (index > index_e)
        return;
    int univ_ij = univ_ij_sort_list(index);

    int ix_natom = univ_ix_natom(univ_ij);
    int iy_natom = univ_iy_natom(univ_ij);
    if (ix_natom * iy_natom <= 0)
        return;

    const int i = univ_cell_pairlist1(1, univ_ij);
    const int j = univ_cell_pairlist1(2, univ_ij);
    const int start_i = start_atom(i);
    const int start_j = start_atom(j);
    int k;

#define sumval(Z) _sumval[(Z) - 1]
    double _sumval[3];
    sumval(1) = 0.0; // virial(1)
    sumval(2) = 0.0; // virial(2)
    sumval(3) = 0.0; // virial(3)
    int8_t check_virial_ij = 0;
    if ((check_virial != 0) && (virial_check(j, i) != 0))
    {
        check_virial_ij = 1;
    }

    if (univ_ij > univ_ncell_near)
    {
        /* far */
        //
        for (int iiy_s = 1 + tile * max_iy_natom; iiy_s <= iy_natom; iiy_s += max_iy_natom << inter_tile_log2)
        { // blocking
            int iiy_e = iiy_s + max_iy_natom - 1;
            if (iiy_e > iy_natom)
                iiy_e = iy_natom;

            // initialize force_iy at shared memory
            for (k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x)
            {
                int n = (k % 3) + 1;
                int iiy = (k / 3) + iiy_s;
                force_iy_smem(n, iiy - iiy_s + 1) = 0.0;
            }

            {
                for (int offset=id_thread_x; offset<iiy_e-iiy_s+1; offset+=num_thread_x) {
                    const int iy=__ldg(&univ_iy_list(iiy_s+offset,univ_ij));
                    const int atom=start_j+iy;
                    j_ids[offset]=iy;
                    j_classes[offset]=__ldg(&atmcls_pbc(atom));
                    for (int d=0; d<4; ++d)
                        j_values[d*max_iy_natom+offset]=__ldg(&coord_pbc(atom,d+1));
                }
                __syncwarp();
            }
            const int iix_max = ix_natom;
            const int iiy_max = iiy_e;
            const int iix_nit = ((iix_max + num_thread_xx - 1) / num_thread_xx) * num_thread_xx;
            const int iiy_nit = ((iiy_max - iiy_s + num_thread_xy) / num_thread_xy) * num_thread_xy + iiy_s - 1;

            for (int iix = id_thread_xx + 1; iix <= iix_nit; iix += num_thread_xx)
            {
                int ix = 0;
                int ixx = 0;
                int iatmcls = 0;
                REAL rtmp1 = 0.0;
                REAL rtmp2 = 0.0;
                REAL rtmp3 = 0.0;
                REAL iqtmp = 0.0;

                if (iix <= iix_max)
                {
                    ix = __ldg(&univ_ix_list(iix, univ_ij));
                    ixx = ix + start_i;
                    iatmcls = __ldg(&atmcls_pbc(ixx));
                    rtmp1 = __ldg(&coord_pbc(ixx, 1)) + __ldg(&cell_move(1, j, i)) * system_x;
                    rtmp2 = __ldg(&coord_pbc(ixx, 2)) + __ldg(&cell_move(2, j, i)) * system_y;
                    rtmp3 = __ldg(&coord_pbc(ixx, 3)) + __ldg(&cell_move(3, j, i)) * system_z;
                    iqtmp = __ldg(&coord_pbc(ixx, 4));
                }

                force_local(1) = 0.0;
                force_local(2) = 0.0;
                force_local(3) = 0.0;

                for (int iiy = id_thread_xy + iiy_s; iiy <= iiy_nit; iiy += num_thread_xy)
                {
                    int iy = 0;
                    int iyy = 0;
                    REAL grad_coef = 0.0;
                    REAL dij1 = 0.0;
                    REAL dij2 = 0.0;
                    REAL dij3 = 0.0;
                    REAL rij2 = 0.0;

                    if (iix <= iix_max && iiy <= iiy_max)
                    {
                        iy = j_ids[iiy-iiy_s];
                        iyy = start_j + iy;

                        dij1 = rtmp1 - j_value(1);
                        dij2 = rtmp2 - j_value(2);
                        dij3 = rtmp3 - j_value(3);
                        rij2 = dij1 * dij1 + dij2 * dij2 + dij3 * dij3;

#if defined(_SINGLE) || defined(_MIXED)
                        rij2=cutoff2*density*__frcp_rn(rij2);
#else
                        rij2=cutoff2*density*__drcp_rn(rij2);
#endif

                        REAL jqtmp = j_value(4);
                        int jatmcls = j_classes[iiy-iiy_s];
                        REAL lj12=__ldg(&nonb_lj12(jatmcls,iatmcls));
                        REAL lj6=__ldg(&nonb_lj6(jatmcls,iatmcls));

                        int L = int(rij2);
                        REAL R = rij2 - L;

                        int L1 = 3 * L - 2;
                        REAL tg0=__ldg(&table_grad(L1)),tg1=__ldg(&table_grad(L1+1));
                        REAL tg2=__ldg(&table_grad(L1+2)),tg3=__ldg(&table_grad(L1+3));
                        REAL tg4=__ldg(&table_grad(L1+4)),tg5=__ldg(&table_grad(L1+5));
                        REAL term_lj12 = tg0 + R * (tg3 - tg0);
                        REAL term_lj6 = tg1 + R * (tg4 - tg1);
                        REAL term_elec = tg2 + R * (tg5 - tg2);

                        grad_coef = term_lj12 * lj12 - term_lj6 * lj6 + iqtmp * jqtmp * term_elec;
                    }

                    REAL work1 = grad_coef * dij1;
                    REAL work2 = grad_coef * dij2;
                    REAL work3 = grad_coef * dij3;

                    force_local(1) -= work1;
                    force_local(2) -= work2;
                    force_local(3) -= work3;

                    if (num_thread_xy == 1) work1 += __shfl_xor_sync(0xffffffff, work1, 1, 32);
                    if (num_thread_xy <= 2) work1 += __shfl_xor_sync(0xffffffff, work1, 2, 32);
                    if (num_thread_xy <= 4) work1 += __shfl_xor_sync(0xffffffff, work1, 4, 32);
                    if (num_thread_xy <= 8) work1 += __shfl_xor_sync(0xffffffff, work1, 8, 32);
                    work1 += __shfl_xor_sync(0xffffffff, work1, 16, 32);

                    if (num_thread_xy == 1) work2 += __shfl_xor_sync(0xffffffff, work2, 1, 32);
                    if (num_thread_xy <= 2) work2 += __shfl_xor_sync(0xffffffff, work2, 2, 32);
                    if (num_thread_xy <= 4) work2 += __shfl_xor_sync(0xffffffff, work2, 4, 32);
                    if (num_thread_xy <= 8) work2 += __shfl_xor_sync(0xffffffff, work2, 8, 32);
                    work2 += __shfl_xor_sync(0xffffffff, work2, 16, 32);

                    if (num_thread_xy == 1) work3 += __shfl_xor_sync(0xffffffff, work3, 1, 32);
                    if (num_thread_xy <= 2) work3 += __shfl_xor_sync(0xffffffff, work3, 2, 32);
                    if (num_thread_xy <= 4) work3 += __shfl_xor_sync(0xffffffff, work3, 4, 32);
                    if (num_thread_xy <= 8) work3 += __shfl_xor_sync(0xffffffff, work3, 8, 32);
                    work3 += __shfl_xor_sync(0xffffffff, work3, 16, 32);

                    if (id_thread_xx == 0 && iix <= iix_max && iiy <= iiy_max)
                    {
                        force_iy_smem(1, iiy - iiy_s + 1) += work1;
                        force_iy_smem(2, iiy - iiy_s + 1) += work2;
                        force_iy_smem(3, iiy - iiy_s + 1) += work3;
                    }
                }

                // update virial
                if (check_virial_ij != 0)
                {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

                if (num_thread_xy >= 2) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 1, 32);
                if (num_thread_xy >= 4) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 2, 32);
                if (num_thread_xy >= 8) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 4, 32);
                if (num_thread_xy >= 16) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 8, 32);

                if (num_thread_xy >= 2) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 1, 32);
                if (num_thread_xy >= 4) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 2, 32);
                if (num_thread_xy >= 8) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 4, 32);
                if (num_thread_xy >= 16) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 8, 32);

                if (num_thread_xy >= 2) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 1, 32);
                if (num_thread_xy >= 4) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 2, 32);
                if (num_thread_xy >= 8) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 4, 32);
                if (num_thread_xy >= 16) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 8, 32);

                if (id_thread_xy == 0 && iix <= iix_max)
                {
                    if (force_local(1) != 0.0)
                        atomicAdd(&(force(ixx, 1)), force_local(1));
                    if (force_local(2) != 0.0)
                        atomicAdd(&(force(ixx, 2)), force_local(2));
                    if (force_local(3) != 0.0)
                        atomicAdd(&(force(ixx, 3)), force_local(3));
                }
            }

            __syncwarp();
            for (k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x)
            {
                int n = (k % 3) + 1;
                int iiy = (k / 3) + iiy_s;
                int iy = __ldg(&univ_iy_list(iiy, univ_ij));
                int iyy = iy + start_j;
                REAL val = force_iy_smem(n, iiy - iiy_s + 1);
                if (val != 0.0)
                    atomicAdd(&(force(iyy, n)), val);
            }
            __syncwarp();
        }
    }
    else
    {
        /* near */
        //
        for (int iiy_s = 1 + tile * max_iy_natom; iiy_s <= iy_natom; iiy_s += max_iy_natom << inter_tile_log2)
        { // blocking
            int iiy_e = iiy_s + max_iy_natom - 1;
            if (iiy_e > iy_natom)
                iiy_e = iy_natom;

            // initialize force_iy at shared memory
            for (k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x)
            {
                int n = (k % 3) + 1;
                int iiy = (k / 3) + iiy_s;
                force_iy_smem(n, iiy - iiy_s + 1) = 0.0;
            }

            {
                for (int offset=id_thread_x; offset<iiy_e-iiy_s+1; offset+=num_thread_x) {
                    const int iy=__ldg(&univ_iy_list(iiy_s+offset,univ_ij));
                    const int atom=start_j+iy;
                    j_ids[offset]=iy;
                    j_classes[offset]=__ldg(&atmcls_pbc(atom));
                    for (int d=0; d<4; ++d)
                        j_values[d*max_iy_natom+offset]=__ldg(&coord_pbc(atom,d+1));
                }
                __syncwarp();
            }
            const int iix_max = ix_natom;
            const int iiy_max = iiy_e;
            const int iix_nit = ((iix_max + num_thread_xx - 1) / num_thread_xx) * num_thread_xx;
            const int iiy_nit = ((iiy_max - iiy_s + num_thread_xy) / num_thread_xy) * num_thread_xy + iiy_s - 1;

            for (int iix = id_thread_xx + 1; iix <= iix_nit; iix += num_thread_xx)
            {
                int ix = 0;
                int ixx = 0;
                int iatmcls = 0;
                REAL rtmp1 = 0.0;
                REAL rtmp2 = 0.0;
                REAL rtmp3 = 0.0;
                REAL iqtmp = 0.0;

                if (iix <= iix_max)
                {
                    ix = __ldg(&univ_ix_list(iix, univ_ij));
                    ixx = ix + start_i;
                    iatmcls = __ldg(&atmcls_pbc(ixx));
                    rtmp1 = __ldg(&coord_pbc(ixx, 1)) + __ldg(&cell_move(1, j, i)) * system_x;
                    rtmp2 = __ldg(&coord_pbc(ixx, 2)) + __ldg(&cell_move(2, j, i)) * system_y;
                    rtmp3 = __ldg(&coord_pbc(ixx, 3)) + __ldg(&cell_move(3, j, i)) * system_z;
                    iqtmp = __ldg(&coord_pbc(ixx, 4));
                }

                force_local(1) = 0.0;
                force_local(2) = 0.0;
                force_local(3) = 0.0;

                int idx0 = (ix - 1) * univ_natom_max;

                for (int iiy = id_thread_xy + iiy_s; iiy <= iiy_nit; iiy += num_thread_xy)
                {
                    int iy = 0;
                    int iyy = 0;
                    REAL grad_coef = 0.0;
                    REAL dij1 = 0.0;
                    REAL dij2 = 0.0;
                    REAL dij3 = 0.0;
                    REAL rij2 = 0.0;

                    if (iix <= iix_max && iiy <= iiy_max)
                    {
                        iy = j_ids[iiy-iiy_s];
                        iyy = start_j + iy;
                        int idx = iy + idx0;
                        // Read the uploaded bit layout directly to keep exclusions compact.
                        const size_t mask_bit = (size_t)(idx - 1) + (size_t)univ_mask2_size * (univ_ij - 1);
                        const uint8_t mask_byte = __ldg((const uint8_t *)_univ_mask2 + (mask_bit >> 3));
                        if ((mask_byte >> (mask_bit & 7)) & 1)
                        {

                            dij1 = rtmp1 - j_value(1);
                            dij2 = rtmp2 - j_value(2);
                            dij3 = rtmp3 - j_value(3);
                            rij2 = dij1 * dij1 + dij2 * dij2 + dij3 * dij3;

#if defined(_SINGLE) || defined(_MIXED)
                            rij2=cutoff2*density*__frcp_rn(rij2);
#else
                            rij2=cutoff2*density*__drcp_rn(rij2);
#endif

                            REAL jqtmp = j_value(4);
                            int jatmcls = j_classes[iiy-iiy_s];
                            REAL lj12=__ldg(&nonb_lj12(jatmcls,iatmcls));
                            REAL lj6=__ldg(&nonb_lj6(jatmcls,iatmcls));

                            int L = int(rij2);
                            REAL R = rij2 - L;

                            int L1 = 3 * L - 2;
                            REAL tg0=__ldg(&table_grad(L1)),tg1=__ldg(&table_grad(L1+1));
                            REAL tg2=__ldg(&table_grad(L1+2)),tg3=__ldg(&table_grad(L1+3));
                            REAL tg4=__ldg(&table_grad(L1+4)),tg5=__ldg(&table_grad(L1+5));
                            REAL term_lj12 = tg0 + R * (tg3 - tg0);
                            REAL term_lj6 = tg1 + R * (tg4 - tg1);
                            REAL term_elec = tg2 + R * (tg5 - tg2);

                            grad_coef = term_lj12 * lj12 - term_lj6 * lj6 + iqtmp * jqtmp * term_elec;
                        }
                    }

                    REAL work1 = grad_coef * dij1;
                    REAL work2 = grad_coef * dij2;
                    REAL work3 = grad_coef * dij3;

                    force_local(1) -= work1;
                    force_local(2) -= work2;
                    force_local(3) -= work3;

                    if (num_thread_xy == 1) work1 += __shfl_xor_sync(0xffffffff, work1, 1, 32);
                    if (num_thread_xy <= 2) work1 += __shfl_xor_sync(0xffffffff, work1, 2, 32);
                    if (num_thread_xy <= 4) work1 += __shfl_xor_sync(0xffffffff, work1, 4, 32);
                    if (num_thread_xy <= 8) work1 += __shfl_xor_sync(0xffffffff, work1, 8, 32);
                    work1 += __shfl_xor_sync(0xffffffff, work1, 16, 32);

                    if (num_thread_xy == 1) work2 += __shfl_xor_sync(0xffffffff, work2, 1, 32);
                    if (num_thread_xy <= 2) work2 += __shfl_xor_sync(0xffffffff, work2, 2, 32);
                    if (num_thread_xy <= 4) work2 += __shfl_xor_sync(0xffffffff, work2, 4, 32);
                    if (num_thread_xy <= 8) work2 += __shfl_xor_sync(0xffffffff, work2, 8, 32);
                    work2 += __shfl_xor_sync(0xffffffff, work2, 16, 32);

                    if (num_thread_xy == 1) work3 += __shfl_xor_sync(0xffffffff, work3, 1, 32);
                    if (num_thread_xy <= 2) work3 += __shfl_xor_sync(0xffffffff, work3, 2, 32);
                    if (num_thread_xy <= 4) work3 += __shfl_xor_sync(0xffffffff, work3, 4, 32);
                    if (num_thread_xy <= 8) work3 += __shfl_xor_sync(0xffffffff, work3, 8, 32);
                    work3 += __shfl_xor_sync(0xffffffff, work3, 16, 32);

                    if (id_thread_xx == 0 && iix <= iix_max && iiy <= iiy_max)
                    {
                        force_iy_smem(1, iiy - iiy_s + 1) += work1;
                        force_iy_smem(2, iiy - iiy_s + 1) += work2;
                        force_iy_smem(3, iiy - iiy_s + 1) += work3;
                    }
                }

                // update virial
                if (check_virial_ij != 0)
                {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

                if (num_thread_xy >= 2) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 1, 32);
                if (num_thread_xy >= 4) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 2, 32);
                if (num_thread_xy >= 8) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 4, 32);
                if (num_thread_xy >= 16) _force_local[0] += __shfl_xor_sync(0xffffffff, _force_local[0], 8, 32);

                if (num_thread_xy >= 2) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 1, 32);
                if (num_thread_xy >= 4) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 2, 32);
                if (num_thread_xy >= 8) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 4, 32);
                if (num_thread_xy >= 16) _force_local[1] += __shfl_xor_sync(0xffffffff, _force_local[1], 8, 32);

                if (num_thread_xy >= 2) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 1, 32);
                if (num_thread_xy >= 4) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 2, 32);
                if (num_thread_xy >= 8) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 4, 32);
                if (num_thread_xy >= 16) _force_local[2] += __shfl_xor_sync(0xffffffff, _force_local[2], 8, 32);

                if (id_thread_xy == 0 && iix <= iix_max)
                {
                    if (force_local(1) != 0.0)
                        atomicAdd(&(force(ixx, 1)), force_local(1));
                    if (force_local(2) != 0.0)
                        atomicAdd(&(force(ixx, 2)), force_local(2));
                    if (force_local(3) != 0.0)
                        atomicAdd(&(force(ixx, 3)), force_local(3));
                }
            }

            __syncwarp();

            // update force(:,:,j)
            for (k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x)
            {
                int n = (k % 3) + 1;
                int iiy = (k / 3) + iiy_s;
                int iy = __ldg(&univ_iy_list(iiy, univ_ij));
                int iyy = iy + start_j;
                REAL val = force_iy_smem(n, iiy - iiy_s + 1);
                if (val != 0.0)
                    atomicAdd(&(force(iyy, n)), val);
            }
            __syncwarp();
        }
    }

    // update virial
    if (check_virial != 0)
    {
        for (int ii = 0; ii < 3; ii++)
        {
            _sumval[ii] += __shfl_xor_sync(0xffffffff, _sumval[ii], 1, 32);
            _sumval[ii] += __shfl_xor_sync(0xffffffff, _sumval[ii], 2, 32);
            _sumval[ii] += __shfl_xor_sync(0xffffffff, _sumval[ii], 4, 32);
            _sumval[ii] += __shfl_xor_sync(0xffffffff, _sumval[ii], 8, 32);
            _sumval[ii] += __shfl_xor_sync(0xffffffff, _sumval[ii], 16, 32);
        }

        // Constant indices keep the three virial accumulators in registers.
        if (id_thread_x == 0) {
            const double value = sumval(1) * (__ldg(&cell_move(1, j, i)) * system_x);
            if (value != 0.0) atomicAdd(&ene_viri_mid(1, index), value);
        }
        if (id_thread_x == 1) {
            const double value = sumval(2) * (__ldg(&cell_move(2, j, i)) * system_y);
            if (value != 0.0) atomicAdd(&ene_viri_mid(2, index), value);
        }
        if (id_thread_x == 2) {
            const double value = sumval(3) * (__ldg(&cell_move(3, j, i)) * system_z);
            if (value != 0.0) atomicAdd(&ene_viri_mid(3, index), value);
        }
    }
}
#undef j_value

namespace genesis {
int nonbond_inter_mode() {
    static const int mode=[] {
        const char *text=std::getenv("GENESIS_NONBOND_INTER");
        if(!text)return 0;
        char *end=nullptr;
        const long value=std::strtol(text,&end,10);
        if(end==text || *end || (value!=0 && value!=38)) {
            std::fprintf(stderr,"GENESIS_NONBOND_INTER must be 0 or 38\n");
            std::abort();
        }
        if(value==38 && sizeof(REAL)!=sizeof(double)) {
            std::fprintf(stderr,"GENESIS_NONBOND_INTER=38 requires double precision\n");
            std::abort();
        }
        return static_cast<int>(value);
    }();
    return mode;
}

void launch_nonbond_inter(dim3 grid,dim3 block,size_t shared_bytes,cudaStream_t stream,int mode,

    const REAL *_coord_pbc,          // ( 1:atom_domain, 1:4 )
    REAL *_force,                    // ( 1:atom_domain, 1:3 )
    double *_ene_virial,             // ( 5 )
    double *_ene_viri_mid,           // ( 1:5, 1:univ_maxcell)
    const int8_t *_cell_move,        // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int *_atmcls_pbc,          // ( 1:atom_domain )
    const int *_natom,               // ( 1:ncel_max )
    const int *_start_atom,          // ( 1:ncel_max )
    const REAL *_nonb_lj12,          // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL *_nonb_lj6,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL *_table_grad,         // ( 1:6*cutoff_int )
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
    REAL density,
    REAL cutoff2,
    REAL pairlistdist2,
    REAL system_x,
    REAL system_y,
    REAL system_z,
    int inter_tile_log2) {
    if(mode!=38)std::abort();
    if(index_e<index_s)return;
    max_iy_natom=80;
    inter_tile_log2=0;
    grid.x=(index_e-index_s+1+block.y-1)/block.y;
    shared_bytes=block.y*max_iy_natom*(7*sizeof(REAL)+2*sizeof(int));
    nonbond_inter_staged_force<<<grid,block,shared_bytes,stream>>>(
        _coord_pbc, _force, _ene_virial, _ene_viri_mid, _cell_move,
        _atmcls_pbc, _natom, _start_atom, _nonb_lj12, _nonb_lj6, _table_grad,
        _univ_cell_pairlist1, _univ_mask2, _univ_ix_natom, _univ_ix_list,
        _univ_iy_natom, _univ_iy_list, _univ_ij_sort_list, _virial_check,
        atom_domain, MaxAtom, MaxAtomCls, num_atom_cls, ncel_local, ncel_bound,
        ncel_max, cutoff_int, univ_maxcell, univ_maxcell1, univ_ncell_near,
        univ_mask2_size, univ_natom_max, index_s, index_e, max_iy_natom,
        check_virial, density, cutoff2, pairlistdist2, system_x, system_y,
        system_z, inter_tile_log2);
    CUDA_CALL(cudaGetLastError());
}
}
