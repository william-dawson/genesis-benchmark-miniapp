#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module improper_outgoing_mod
 use bonded_incoming_mod,only:bonded_improper_bootstrap
 use iso_c_binding,only:c_ptr,c_int,c_loc
 use sp_domain_str_mod
 use sp_enefunc_str_mod
 use mpi_parallel_mod
 implicit none
 private
 public::improper_outgoing,improper_outgoing_check,improper_outgoing_verify
 interface
  integer(c_int) function mode_c() bind(C,name='genesis_improper_outgoing_mode')
   import c_int
  end function
  integer(c_int) function run_c(owner,pitch,total,verify) bind(C,name='genesis_improper_outgoing')
   import c_ptr,c_int
   type(c_ptr),value::owner
   integer(c_int),value::pitch,total,verify
  end function
  subroutine check_c() bind(C,name='genesis_improper_outgoing_check')
  end subroutine
 end interface
contains
 logical function improper_outgoing_verify()
  improper_outgoing_verify=mode_c()==2
 end function
 logical function improper_outgoing(d,e)
  type(s_domain),target,intent(in)::d
  type(s_enefunc),target,intent(inout)::e
  improper_outgoing=.false.
#ifdef USE_GPU
  if(mode_c()==0)return
  if(nproc_country/=1.or.d%num_cell_boundary/=0.or.e%forcefield/=ForcefieldCHARMM) &
    error stop 'GPU improper outgoing requires one-rank CHARMM without halos'
  call bonded_improper_bootstrap(d,e)
  improper_outgoing=run_c(c_loc(d%cell_pair),size(d%cell_pair,1),e%num_impr_all,merge(1,0,improper_outgoing_verify()))/=0
#endif
 end function
 subroutine improper_outgoing_check()
  call check_c()
 end subroutine
end module
