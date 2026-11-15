#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module maintenance_state_mod
  use iso_c_binding,only:c_int,c_double
  use cudafor
  use sp_domain_str_mod
  use constants_mod
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_nan
  implicit none
  private
  public::maintenance_state_begin,maintenance_state_consume,maintenance_state_poison
  interface
    integer(c_int) function mode_c() bind(C,name='genesis_maint_state_mode')
      import c_int
    end function
    integer(c_int) function begin_c(values,width,cells,eligible) bind(C,name='genesis_maint_state_begin')
      import c_int,c_double
      real(c_double),device,intent(in)::values(*)
      integer(c_int),value::width,cells,eligible
    end function
    integer(c_int) function consume_c(values,first,second,width,cells,short,work,long) &
        bind(C,name='genesis_maint_state_consume')
      import c_int,c_double
      real(c_double),device,intent(inout)::values(*),first(*),second(*)
      integer(c_int),value::width,cells,short,work,long
    end function
  end interface
contains
  logical function maintenance_state_begin(values,width,cells,eligible)
    integer,intent(in)::width,cells
    logical,intent(in)::eligible
    real(8),device,intent(in)::values(:,:,:)
    maintenance_state_begin=begin_c(values,width,cells,merge(1,0,eligible))/=0
  end function

  logical function maintenance_state_consume(domain,values,first,second,width,cells,short,work,long)
    type(s_domain),intent(in)::domain
    integer::cell
    integer,intent(in)::width,cells
    logical,intent(in)::short,work,long
    real(8),device,intent(inout)::values(:,:,:),first(:,:,:),second(:,:,:)
    maintenance_state_consume=consume_c(values,first,second,width,cells, &
      merge(1,0,short),merge(1,0,work),merge(1,0,long))/=0
    if(maintenance_state_consume.and.mode_c()==2) then
      do cell=1,domain%num_cell_local
        if(.not.all(ieee_is_nan(domain%coord_ref(:,1:domain%num_atom(cell),cell)))) &
          error stop 'Retired CPU coordinate reference was overwritten'
        if(.not.all(ieee_is_nan(domain%velocity_full(:,1:domain%num_atom(cell),cell)))) &
          error stop 'Retired CPU work was overwritten'
        if(.not.all(ieee_is_nan(domain%velocity_ref(:,1:domain%num_atom(cell),cell)))) &
          error stop 'Retired CPU velocity reference was overwritten'
      end do
      write(*,'(a)') 'GENESIS_MAINT_STATE: CPU reference mirrors remain poisoned'
    end if
  end function

  subroutine maintenance_state_poison(domain)
    type(s_domain),intent(inout)::domain
    if(mode_c()/=2) return
    domain%coord_ref=ieee_value(0.0_wip,ieee_quiet_nan)
    domain%velocity_full=ieee_value(0.0_wip,ieee_quiet_nan)
    domain%velocity_ref=ieee_value(0.0_wip,ieee_quiet_nan)
    domain%force_short=ieee_value(0.0_wip,ieee_quiet_nan)
    domain%force_long=ieee_value(0.0_wip,ieee_quiet_nan)
    write(*,'(a)') 'GENESIS_MAINT_STATE: five CPU mirrors poisoned'
  end subroutine
end module
