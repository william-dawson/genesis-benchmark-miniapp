module genesis_pair_sort_mod
  use iso_c_binding,only:c_int
  implicit none
  private
  public::pair_sort_enabled,pair_sort_verify,pair_sort_build,pair_sort_check
  interface
    subroutine pair_sort_build(total,local,near,width,maxload,nonzero,ix,iy,loads,order) &
        bind(C,name='genesis_pair_sort_build')
      import c_int
      integer(c_int),value::total,local,near,width,maxload
      integer(c_int),intent(out)::nonzero
      integer(c_int),intent(inout)::ix(*),iy(*),loads(*),order(*)
    end subroutine
    subroutine pair_sort_check(total,nonzero,loads,order) bind(C,name='genesis_pair_sort_check')
      import c_int
      integer(c_int),value::total,nonzero
      integer(c_int),intent(in)::loads(*),order(*)
    end subroutine
  end interface
contains
  logical function flag(name)
    character(*),intent(in)::name
    character(16)::value
    integer::status
    call get_environment_variable(name,value,status=status)
    flag=.false.
    if(status==1)return
    if(status/=0)error stop 'Truncated pair sorting flag'
    if(trim(value)/='0'.and.trim(value)/='1')error stop 'Invalid pair sorting flag'
    flag=trim(value)=='1'
  end function
  logical function pair_sort_enabled()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_PAIR_SORT');initialized=.true.
    end if
    pair_sort_enabled=selected
  end function
  logical function pair_sort_verify()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_PAIR_SORT_VERIFY');initialized=.true.
    end if
    pair_sort_verify=selected
  end function
end module
