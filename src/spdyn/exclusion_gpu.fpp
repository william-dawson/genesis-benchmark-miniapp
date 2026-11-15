#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module exclusion_gpu_mod
 use constants_mod
 use iso_c_binding,only:c_ptr,c_int,c_double,c_loc
 use sp_domain_str_mod
 use sp_enefunc_str_mod
 use sp_constraints_str_mod
 use mpi_parallel_mod
 use genesis_domain_masks_mod,only:domain_mask_mode
 implicit none
 private
 public::exclusion_gpu
 type,bind(C)::view
  type(c_ptr)::p(27)
  integer(c_int)::n(20)
 end type
 interface
  integer(c_int) function enabled_c() bind(C,name='genesis_exclusion_enabled')
   import c_int
  end function
  subroutine build_c(v,epsilon,maximum,total) bind(C,name='genesis_exclusion_build')
   import view,c_double,c_int
   type(view),intent(in)::v
   real(c_double),value::epsilon
   integer(c_int),intent(out)::maximum,total
  end subroutine
 end interface
contains
 logical function exclusion_gpu(first,constraint,c,d,e,reference)
 logical,intent(in)::first,constraint,reference
 type(s_constraints),target,intent(in)::c
 type(s_domain),target,intent(inout)::d
 type(s_enefunc),target,intent(inout)::e
 type(view)::v
 exclusion_gpu=.false.
#ifdef USE_GPU
 if(first.or.reference.or.enabled_c()==0)return
 if(nproc_country/=1.or.d%num_cell_boundary/=0.or.storage_size(1.0_wp)/=64.or. &
    e%forcefield/=ForcefieldCHARMM.or.d%nonbond_kernel/=NBK_GPU.or.domain_mask_mode()/=6) &
    error stop 'GPU exclusion requires one-rank CHARMM double GPU masks mode6'
 if(c%water_type==TIP4.or.sum(e%num_rb_dihedral)>0)error stop 'GPU exclusion requires TIP3 and no RB dihedrals'
 v%p(1)=c_loc(d%num_atom)
 v%p(2)=c_loc(d%num_water)
 v%p(3)=c_loc(d%water_list)
 v%p(4)=c_loc(d%id_g2l)
 v%p(5)=c_loc(d%cell_pair)
 v%p(6)=c_loc(d%cell_pairlist2)
 v%p(7)=c_loc(d%charge)
 v%p(8)=c_loc(e%num_bond)
 v%p(9)=c_loc(e%bond_list)
 v%p(10)=c_loc(e%bond_kind)
 v%p(11)=c_loc(e%num_angle)
 v%p(12)=c_loc(e%angle_list)
 v%p(13)=c_loc(e%angle_kind)
 v%p(14)=c_loc(e%num_dihedral)
 v%p(15)=c_loc(e%dihe_list)
 v%p(16)=c_loc(e%dihe_kind)
 v%p(17)=c_loc(e%dihe_periodicity)
 v%p(18)=c_loc(c%HGr_local)
 v%p(19)=c_loc(c%HGr_bond_list)
 v%p(20)=c_loc(e%num_nonb_excl)
 v%p(21)=c_loc(e%num_nb14_calc)
 v%p(22)=c_loc(e%nonb_excl_list)
 v%p(23)=c_loc(e%nb14_calc_list)
 v%p(24)=c_loc(e%sc_calc_list)
 v%p(25)=c_loc(d%start_atom)
 v%p(26)=c_loc(e%exclusion_mask1)
 v%p(27)=c_loc(e%exclusion_mask)
 v%n(1)=d%num_cell_local
 v%n(2)=MaxAtom
 v%n(3)=size(e%exclusion_mask,3)
 v%n(4)=size(d%water_list,1)
 v%n(5)=size(d%water_list,2)
 v%n(6)=size(d%cell_pair,1)
 v%n(7)=size(d%cell_pairlist2,1)
 v%n(8)=size(d%charge,1)
 v%n(9)=size(e%bond_list,2)
 v%n(10)=size(e%angle_list,2)
 v%n(11)=size(e%dihe_list,2)
 v%n(12)=c%connect
 v%n(13)=size(c%HGr_bond_list,1)
 v%n(14)=size(c%HGr_bond_list,2)
 v%n(15)=size(e%nonb_excl_list,2)
 v%n(16)=size(e%nb14_calc_list,2)
 v%n(17)=e%notation_14types
 v%n(18)=e%excl_level
 v%n(19)=merge(1,0,constraint)
 v%n(20)=size(d%id_g2l,2)
 call build_c(v,EPS,d%max_num_atom,d%num_atom_domain)
 exclusion_gpu=.true.
#endif
 end function
end module
