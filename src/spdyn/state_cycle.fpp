#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module state_cycle_mod
 use iso_c_binding,only:c_int,c_double
 use constants_mod,only:wip
 use cudafor
 use sp_domain_str_mod
 implicit none
 private
 public::state_cycle_begin,state_cycle_consume,state_cycle_retired
 interface
  integer(c_int) function begin_c(x,v,hx,hv,n,width,cells,eligible) bind(C,name='genesis_state_cycle_begin')
   import c_int,wip
   real(wip),device,intent(in)::x(*),v(*)
   real(wip),intent(inout)::hx(*),hv(*)
   integer(c_int),intent(in)::n(*)
   integer(c_int),value::width,cells,eligible
  end function
  integer(c_int) function consume_c(x,v,hx,hv,n,width,cells,eligible) bind(C,name='genesis_state_cycle_consume')
   import c_int,wip
   real(wip),device,intent(in)::x(*),v(*)
   real(wip),intent(in)::hx(*),hv(*)
   integer(c_int),intent(in)::n(*)
   integer(c_int),value::width,cells,eligible
  end function
  integer(c_int) function state_cycle_retired() bind(C,name='genesis_state_cycle_retired')
   import c_int
  end function
 end interface
contains
 logical function state_cycle_begin(d,x,v,eligible) result(skip)
  type(s_domain),intent(inout)::d
  real(wip),device,intent(in)::x(:,:,:),v(:,:,:)
  logical,intent(in)::eligible
  skip=.false.
  skip=begin_c(x,v,d%coord,d%velocity,d%num_atom,size(x,2),size(x,3),merge(1,0,eligible))/=0
 end function
 logical function state_cycle_consume(d,x,v,eligible) result(retained)
  type(s_domain),intent(in)::d
  real(wip),device,intent(in)::x(:,:,:),v(:,:,:)
  logical,intent(in)::eligible
  retained=.false.
  retained=consume_c(x,v,d%coord,d%velocity,d%num_atom,size(x,2),size(x,3),merge(1,0,eligible))/=0
 end function
end module
