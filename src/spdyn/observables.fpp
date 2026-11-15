#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module observables_mod
 use iso_c_binding,only:c_int,c_double
 use cudafor
 use sp_domain_str_mod
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_nan
 implicit none
 private
 public::observables_stage,observables_take,observables_complete,observables_poison,observables_mode
 logical,save::poison_pending=.false.
 interface
  integer(c_int) function observables_mode() bind(C,name='genesis_observables_mode')
   import c_int
  end function
  integer(c_int) function stage_c(f,v,m,n,width,cells,step) bind(C,name='genesis_observables_stage')
   import c_int,c_double
   real(c_double),device,intent(in)::f(*),v(*),m(*)
   integer(c_int),device,intent(in)::n(*)
   integer(c_int),value::width,cells,step
  end function
  integer(c_int) function take_c(step,rmsg,kinetic) bind(C,name='genesis_observables_take')
   import c_int,c_double
   integer(c_int),value::step
   real(c_double),intent(out)::rmsg,kinetic
  end function
  subroutine check_c(step,rmsg,kinetic) bind(C,name='genesis_observables_check')
   import c_int,c_double
   integer(c_int),value::step
   real(c_double),value::rmsg,kinetic
  end subroutine
 end interface
contains
 logical function observables_stage(f,v,m,n,step) result(skip)
  real(8),device,intent(in)::f(:,:,:),v(:,:,:),m(:,:)
  integer,device,intent(in)::n(:)
  integer,intent(in)::step
  skip=stage_c(f,v,m,n,size(f,2),size(f,3),step)==1
 end function
 logical function observables_take(step,rmsg,kinetic) result(taken)
  integer,intent(in)::step
  real(8),intent(out)::rmsg,kinetic
  taken=take_c(step,rmsg,kinetic)/=0
 end function
 subroutine observables_poison(d)
  type(s_domain),intent(inout)::d
  if(observables_mode()/=3)return
  d%coord=ieee_value(0.0d0,ieee_quiet_nan);d%velocity=ieee_value(0.0d0,ieee_quiet_nan)
  d%coord_ref=ieee_value(0.0d0,ieee_quiet_nan);d%velocity_ref=ieee_value(0.0d0,ieee_quiet_nan)
  d%velocity_full=ieee_value(0.0d0,ieee_quiet_nan);d%force_short=ieee_value(0.0d0,ieee_quiet_nan)
  d%force_long=ieee_value(0.0d0,ieee_quiet_nan);d%force=ieee_value(0.0d0,ieee_quiet_nan)
  d%velocity_half=ieee_value(0.0d0,ieee_quiet_nan)
  poison_pending=.true.
  write(*,'(a)')'GENESIS_OBSERVABLES: CPU reporting mirrors poisoned'
 end subroutine
 subroutine observables_complete(d,step,rmsg,kinetic)
  type(s_domain),intent(in)::d
  integer,intent(in)::step
  real(8),intent(in)::rmsg,kinetic
  integer::c,n
  call check_c(step,rmsg,kinetic)
  if(.not.poison_pending)return
  do c=1,d%num_cell_local
   n=d%num_atom(c)
   if(.not.all(ieee_is_nan(d%force(:,1:n,c))).or..not.all(ieee_is_nan(d%velocity_ref(:,1:n,c)))) &
    error stop 'Retired reporting CPU inputs were rewritten'
  end do
  poison_pending=.false.
  write(*,'(a)')'GENESIS_OBSERVABLES: scalar report and poisoned CPU inputs pass'
 end subroutine
end module
