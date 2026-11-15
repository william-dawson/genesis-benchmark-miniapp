#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module final_reductions_mod
  use iso_c_binding,only:c_int,c_double,c_int64_t
  use cudafor
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  public::final_reductions_enabled,final_constraint_reduce,final_thermostat_reduce
  interface
    integer(c_int) function tree_c() bind(C,name="genesis_final_reductions_tree")
      import c_int
    end function
    integer(c_int) function mode_c() bind(C,name="genesis_final_reductions_mode")
      import c_int
    end function
    subroutine reduce_checked_c(input,output,cells,failure,bad) bind(C,name="genesis_final_reduce_checked")
      import c_int,c_double
      real(c_double),device,intent(in)::input(*)
      real(c_double),intent(out)::output(*)
      integer(c_int),value::cells
      integer(c_int),device,intent(in)::failure
      integer(c_int),intent(out)::bad
    end subroutine
    subroutine reduce_c(input,output,cells) bind(C,name="genesis_final_reduce")
      import c_int,c_double
      real(c_double),device,intent(in)::input(*)
      real(c_double),intent(out)::output(*)
      integer(c_int),value::cells
    end subroutine
  end interface
contains
  logical function final_reductions_enabled()
    final_reductions_enabled=mode_c()/=0
  end function
  subroutine final_constraint_reduce(input,mirror,output,failure)
    real(8),device,intent(in)::input(:,:,:)
    real(8),intent(inout)::mirror(:,:,:)
    real(8),intent(out)::output(3,3)
    integer(c_int),device,intent(in)::failure
    real(8)::expected(3,3)
    integer::cell,bad
    if(any(shape(input)/=shape(mirror)).or.size(input,1)/=3.or.size(input,2)/=3) &
      error stop 'Invalid constraint final reduction shape'
    call reduce_checked_c(input,output,size(input,3),failure,bad)
    if(bad/=0)error stop 'GPU SHAKE failed to converge'
    if(.not.all(ieee_is_finite(output)))error stop 'Nonfinite constraint total'
    if(mode_c()/=2)return
    mirror=input;expected=0
    do cell=1,size(input,3)
      expected=expected+mirror(:,:,cell)
    end do
    if(tree_c()/=0)then
      call check_budget(reshape(output,[9]),reshape(expected,[9]),1.0d-7,'constraint')
    else
      if(any(transfer(output,[0_c_int64_t],9)/=transfer(expected,[0_c_int64_t],9))) &
        error stop 'Constraint final sum changed bits'
      write(*,'(a)') 'GENESIS_FINAL_REDUCTIONS: exact constraint total and mirror poison pass'
    end if
    mirror=ieee_value(0.0d0,ieee_quiet_nan)
  end subroutine
  subroutine final_thermostat_reduce(input,mirror,output)
    real(8),device,intent(in)::input(:,:)
    real(8),intent(inout)::mirror(:,:)
    real(8),intent(out)::output(9)
    real(8)::expected(9)
    integer::cell
    if(any(shape(input)/=shape(mirror)).or.size(input,1)/=9)error stop 'Invalid thermostat final reduction shape'
    call reduce_c(input,output,size(input,2))
    if(.not.all(ieee_is_finite(output)))error stop 'Nonfinite thermostat total'
    if(mode_c()/=2)return
    mirror=input;expected=0
    do cell=1,size(input,2)
      expected=expected+mirror(:,cell)
    end do
    if(tree_c()/=0)then
      call check_budget(output,expected,1.0d-8,'thermostat')
    else
      if(any(transfer(output,[0_c_int64_t],9)/=transfer(expected,[0_c_int64_t],9))) &
        error stop 'Thermostat final sum changed bits'
      write(*,'(a)') 'GENESIS_FINAL_REDUCTIONS: exact thermostat total and mirror poison pass'
    end if
    mirror=ieee_value(0.0d0,ieee_quiet_nan)
  end subroutine
  subroutine check_budget(actual,expected,atol,label)
    real(8),intent(in)::actual(9),expected(9),atol
    character(*),intent(in)::label
    real(8)::ratio,rms
    if(.not.all(ieee_is_finite(expected)))error stop 'Nonfinite final sum reference'
    ratio=maxval(abs(actual-expected)/(atol+1.0d-10*abs(expected)))
    rms=sqrt(sum((actual-expected)**2))/max(sqrt(sum(expected**2)),tiny(1.0d0))
    if(ratio>1.or.rms>1.0d-11)error stop 'GPU final sum exceeds unchanged numerical budget'
    write(*,'(a,a,a,2(es20.12,1x))') 'GENESIS_FINAL_REDUCTIONS: budgeted ',label,' total and mirror poison pass ',ratio,rms
  end subroutine
end module
