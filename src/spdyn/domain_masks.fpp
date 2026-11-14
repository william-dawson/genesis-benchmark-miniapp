module genesis_domain_masks_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_loc, c_int, c_size_t, c_int8_t, c_int16_t
  implicit none
  private
  public :: domain_mask_mode, domain_fill_masks
  interface
    function host_memset(address, value, bytes) bind(c,name='memset') result(destination)
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: address
      integer(c_int), value :: value
      integer(c_size_t), value :: bytes
      type(c_ptr) :: destination
    end function
  end interface
contains
  integer function domain_mask_mode()
    character(16) :: value
    integer :: status, code
    integer, save :: selected = -1
    if (selected < 0) then
      selected = 0
      call get_environment_variable('GENESIS_DOMAIN_MASK', value, status=status)
      if (status == 0) then
        read(value, *, iostat=code) selected
        if (code /= 0) error stop 'Invalid GENESIS_DOMAIN_MASK'
        if (selected < 0 .or. selected > 6) error stop 'GENESIS_DOMAIN_MASK must be between 0 and 6'
      end if
    end if
    domain_mask_mode = selected
  end function

  subroutine domain_fill_masks(ld, nlocal, npairs, natom, pairs, mask1_address, mask_address, mode)
    integer, intent(in) :: ld, nlocal, npairs, mode
    integer, intent(in) :: natom(*)
    integer(c_int16_t), intent(in) :: pairs(2,*)
    type(c_ptr), value :: mask1_address, mask_address
    integer(c_int8_t), pointer, contiguous :: mask1(:,:,:), mask(:,:,:)
    integer(c_int8_t), pointer, contiguous :: flat1(:), flat(:)
    integer :: cell, pair, row, column, count
    integer :: status
    character(16) :: value
    logical :: verify
    type(c_ptr) :: destination

    call c_f_pointer(mask1_address, mask1, [ld,ld,nlocal])
    call c_f_pointer(mask_address, mask, [ld,ld,npairs])
    call c_f_pointer(mask1_address, flat1, [int(ld,8)*ld*nlocal])
    call c_f_pointer(mask_address, flat, [int(ld,8)*ld*npairs])
    !$omp parallel private(cell,pair,row,column,count,destination)
    !$omp do schedule(static)
    do cell = 1, nlocal
      count = natom(cell)
      if (mode == 2) count = ld
      if (mode == 3) then
        flat1((int(cell,8)-1)*ld*ld+1:int(cell,8)*ld*ld) = 1_c_int8_t
      else if (mode == 4) then
        ! Fill one contiguous span covering the active square, including row padding.
        if (count > 0) flat1((int(cell,8)-1)*ld*ld+1:(int(cell,8)-1)*ld*ld+(count-1)*ld+count) = 1_c_int8_t
      else if (mode == 5) then
        destination = host_memset(c_loc(mask1(1,1,cell)),1_c_int,int(ld,c_size_t)*ld)
      else if (mode == 6) then
        if (count > 0) destination = host_memset(c_loc(mask1(1,1,cell)),1_c_int, &
                                                 int(count-1,c_size_t)*ld+count)
      else
      do column = 1, count
        !$omp simd
        do row = 1, count
          mask1(row,column,cell) = 1_c_int8_t
        end do
      end do
      end if
      do row = 1, natom(cell)
        mask1(row,row,cell) = 0_c_int8_t
      end do
    end do
    !$omp end do
    !$omp do schedule(static)
    do pair = 1, npairs
      count = max(natom(pairs(1,pair)),natom(pairs(2,pair)))
      if (mode == 2) count = ld
      if (mode == 3) then
        flat((int(pair,8)-1)*ld*ld+1:int(pair,8)*ld*ld) = 1_c_int8_t
      else if (mode == 4) then
        if (count > 0) flat((int(pair,8)-1)*ld*ld+1:(int(pair,8)-1)*ld*ld+(count-1)*ld+count) = 1_c_int8_t
      else if (mode == 5) then
        destination = host_memset(c_loc(mask(1,1,pair)),1_c_int,int(ld,c_size_t)*ld)
      else if (mode == 6) then
        if (count > 0) destination = host_memset(c_loc(mask(1,1,pair)),1_c_int, &
                                                 int(count-1,c_size_t)*ld+count)
      else
      do column = 1, count
        !$omp simd
        do row = 1, count
          mask(row,column,pair) = 1_c_int8_t
        end do
      end do
      end if
    end do
    !$omp end do
    !$omp end parallel

    call get_environment_variable('GENESIS_DOMAIN_VERIFY', value, status=status)
    verify = status == 0 .and. trim(value) == '1'
    if (.not. verify) return
    do cell = 1, nlocal
      do column = 1, natom(cell)
        do row = 1, natom(cell)
          if (row == column) then
            if (mask1(row,column,cell) /= 0) error stop 'Domain intra mask diagonal mismatch'
          else
            if (mask1(row,column,cell) /= 1) error stop 'Domain intra mask mismatch'
          end if
        end do
      end do
    end do
    do pair = 1, npairs
      count = max(natom(pairs(1,pair)),natom(pairs(2,pair)))
      if (any(mask(1:count,1:count,pair) /= 1)) error stop 'Domain inter mask mismatch'
    end do
    write(*,'(a,3i12)') 'GENESIS_DOMAIN_VERIFY: mask initialization passes ', nlocal, npairs, ld
  end subroutine
end module
