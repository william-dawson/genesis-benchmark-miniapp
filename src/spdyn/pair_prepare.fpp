#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module pair_prepare_mod
 use constants_mod
 use iso_c_binding,only:c_int,c_ptr,c_loc
 use sp_domain_str_mod
 use mpi_parallel_mod
 implicit none
 private
 public::pair_prepare_begin,pair_prepare_reference
 interface
  integer(c_int) function mode_c() bind(C,name='genesis_pair_prepare_mode')
   import c_int
  end function
  integer(c_int) function begin_c(coords,classes,cells,width,atoms) bind(C,name='genesis_pair_prepare_begin')
   import c_int,c_ptr
   type(c_ptr),value::coords,classes
   integer(c_int),value::cells,width,atoms
  end function
 end interface
contains
 logical function pair_prepare_reference()
  pair_prepare_reference=mode_c()==2
 end function
 logical function pair_prepare_begin(d,maximum) result(active)
  type(s_domain),target,intent(in)::d
  integer,intent(out)::maximum
  active=.false.;maximum=0
#ifdef USE_GPU
  if(mode_c()==0)return
  if(nproc_country/=1.or.d%num_cell_boundary/=0.or.storage_size(1.0_wp)/=64) &
   error stop 'GPU pair preparation requires one-rank double without halos'
  active=begin_c(c_loc(d%translated),c_loc(d%atmcls_pbc),d%num_cell_local,size(d%coord,2),d%num_atom_domain)/=0
  maximum=d%max_num_atom
#endif
 end function
end module
