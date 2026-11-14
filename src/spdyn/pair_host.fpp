module genesis_pair_host_mod
  use iso_c_binding, only: c_ptr, c_f_pointer, c_loc, c_int, c_size_t, c_int8_t
#ifdef OMP
  use omp_lib, only: omp_get_wtime
#endif
  implicit none
  private
  public :: pair_mode, pair_clock, pair_report, pair_fill, pair_pack, pair_verify_fill, pair_verify_pack
  interface
    subroutine copy_rectangle(input, output, rows, columns, input_ld, output_ld) bind(c,name='genesis_pair_copy')
      import :: c_ptr, c_int
      type(c_ptr), value :: input, output
      integer(c_int), value :: rows, columns, input_ld, output_ld
    end subroutine
    function copy_bytes(destination, source, bytes) bind(c,name='memcpy') result(address)
      import :: c_ptr, c_size_t
      type(c_ptr), value :: destination, source
      integer(c_size_t), value :: bytes
      type(c_ptr) :: address
    end function
    subroutine pack_bytes(input, output, count) bind(c,name='genesis_pair_pack')
      import :: c_ptr, c_int
      type(c_ptr), value :: input, output
      integer(c_int), value :: count
    end subroutine
  end interface
contains
  logical function verify_enabled()
    logical, save :: initialized=.false., enabled=.false.
    character(16) :: value
    integer :: status
    if (.not. initialized) then
      call get_environment_variable('GENESIS_PAIR_VERIFY',value,status=status)
      enabled=status==0 .and. trim(value)=='1'
      initialized=.true.
    end if
    verify_enabled=enabled
  end function

  integer function pair_mode()
    integer, save :: selected = -1
    integer :: status, code
    character(16) :: value
    if (selected < 0) then
      selected=0
      call get_environment_variable('GENESIS_PAIR_HOST',value,status=status)
      if (status == 0) then
        read(value,*,iostat=code) selected
        if (code /= 0) error stop 'Invalid GENESIS_PAIR_HOST'
        if (selected < 0 .or. selected > 7) error stop 'GENESIS_PAIR_HOST must be between 0 and 7'
      end if
    end if
    pair_mode=selected
  end function

  real(8) function pair_clock()
    logical, save :: initialized=.false., enabled=.false.
    character(16) :: value
    integer :: status
    integer(8) :: ticks, rate
    if (.not. initialized) then
      call get_environment_variable('GENESIS_PAIR_TIMING',value,status=status)
      enabled=status==0 .and. trim(value)=='1'
      initialized=.true.
    end if
    pair_clock=0.0d0
    if (.not. enabled) return
#ifdef OMP
    pair_clock=omp_get_wtime()
#else
    call system_clock(ticks,rate)
    pair_clock=real(ticks,8)/real(rate,8)
#endif
  end function

  subroutine pair_report(times)
    real(8), intent(in) :: times(6)
    if (times(1)==0.0d0) return
    write(*,'(a,5(es18.10,1x))') 'GENESIS_PAIR_PHASES: ', (times(2:6)-times(1:5))*1.0d6
  end subroutine

  subroutine pair_fill(ld, pitch, width, nlocal, nnear, natom, pairs, one_address, inter_address, out_address)
    integer, intent(in) :: ld, pitch, width, nlocal, nnear, natom(*), pairs(2,*)
    type(c_ptr), value :: one_address, inter_address, out_address
    integer(c_int8_t), pointer, contiguous :: one(:,:,:), inter(:,:,:), output(:,:)
    type(c_ptr) :: address
    integer :: plane, cell, first, second, column, row, count
    call c_f_pointer(one_address,one,[ld,ld,nlocal])
    call c_f_pointer(inter_address,inter,[ld,ld,nnear-nlocal])
    call c_f_pointer(out_address,output,[pitch,nnear])
    !$omp parallel do schedule(static) private(cell,column,row,count,address)
    do plane=1,nlocal
      cell=pairs(1,plane)
      count=natom(cell)
      do column=1,count
        do row=1,column-1
          output((column-1)*width+row,plane)=one(column,row,cell)
        end do
        address=copy_bytes(c_loc(output((column-1)*width+column,plane)), &
                           c_loc(one(column,column,cell)),int(count-column+1,c_size_t))
      end do
    end do
    !$omp end parallel do
    !$omp parallel do schedule(static) private(first,second,column,count,address)
    do plane=nlocal+1,nnear
      first=pairs(1,plane)
      second=pairs(2,plane)
      count=natom(second)
      if (btest(pair_mode(),2)) then
        if (count>0 .and. natom(first)>0) &
          call copy_rectangle(c_loc(inter(1,1,plane-nlocal)),c_loc(output(1,plane)),count,natom(first),ld,width)
      else
        do column=1,natom(first)
          address=copy_bytes(c_loc(output((column-1)*width+1,plane)), &
                             c_loc(inter(1,column,plane-nlocal)),int(count,c_size_t))
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine

  subroutine pair_pack(input_address, output_address, start_index, end_index, output_size)
    type(c_ptr), value :: input_address, output_address
    integer, intent(in) :: start_index, end_index, output_size
    integer(c_int8_t), pointer, contiguous :: input(:), output(:)
    integer :: groups, offset, block, first, count, remainder, bit
    integer, parameter :: chunk=16384
    integer(c_int8_t) :: last
    call c_f_pointer(input_address,input,[end_index])
    call c_f_pointer(output_address,output,[output_size])
    groups=(end_index-start_index+1)/8
    offset=start_index/8
    !$omp parallel do schedule(static) private(first,count)
    do block=0,(groups+chunk-1)/chunk-1
      first=block*chunk
      count=min(chunk,groups-first)
      call pack_bytes(c_loc(input(start_index+8*first)),c_loc(output(offset+first+1)),count)
    end do
    !$omp end parallel do
    remainder=mod(end_index-start_index+1,8)
    if (remainder > 0) then
      last=0
      do bit=0,remainder-1
        if (input(start_index+groups*8+bit)==1) last=ibset(last,bit)
      end do
      output(offset+groups+1)=last
    end if
  end subroutine

  subroutine pair_verify_fill(ld, pitch, width, nlocal, nnear, natom, pairs, one_address, inter_address, out_address)
    integer, intent(in) :: ld, pitch, width, nlocal, nnear, natom(*), pairs(2,*)
    type(c_ptr), value :: one_address, inter_address, out_address
    integer(c_int8_t), pointer, contiguous :: one(:,:,:), inter(:,:,:), output(:,:)
    integer :: plane, first, second, column, row
    integer(c_int8_t) :: expected
    integer(8) :: failures
    if (.not. verify_enabled()) return
    call c_f_pointer(one_address,one,[ld,ld,nlocal])
    call c_f_pointer(inter_address,inter,[ld,ld,nnear-nlocal])
    call c_f_pointer(out_address,output,[pitch,nnear])
    failures=0
    !$omp parallel do schedule(static) private(first,second,column,row,expected) reduction(+:failures)
    do plane=1,nnear
      first=pairs(1,plane)
      second=pairs(2,plane)
      do column=1,natom(first)
        do row=1,natom(second)
          if (plane <= nlocal) then
            expected=one(max(row,column),min(row,column),first)
          else
            expected=inter(row,column,plane-nlocal)
          end if
          if (output((column-1)*width+row,plane)/=expected) failures=failures+1
        end do
      end do
    end do
    !$omp end parallel do
    if (failures/=0) error stop 'Pair byte-mask verification failed'
  end subroutine

  subroutine pair_verify_pack(input_address, output_address, start_index, end_index, output_size, complete)
    type(c_ptr), value :: input_address, output_address
    integer, intent(in) :: start_index, end_index, output_size
    logical, intent(in) :: complete
    integer(c_int8_t), pointer, contiguous :: input(:), output(:)
    integer :: index, packed, bit
    integer(8) :: failures
    if (.not. verify_enabled()) return
    call c_f_pointer(input_address,input,[end_index])
    call c_f_pointer(output_address,output,[output_size])
    failures=0
    !$omp parallel do schedule(static) private(packed,bit) reduction(+:failures)
    do index=start_index,end_index
      packed=(index-1)/8+1
      bit=mod(index-1,8)
      if (btest(output(packed),bit).neqv.(input(index)==1)) failures=failures+1
    end do
    !$omp end parallel do
    if (mod(end_index,8)>0) then
      do bit=mod(end_index,8),7
        if (btest(output((end_index-1)/8+1),bit)) failures=failures+1
      end do
    end if
    if (failures/=0) error stop 'Pair packing verification failed'
    if (complete) write(*,'(a)') 'GENESIS_PAIR_VERIFY: mask and packing pass'
  end subroutine
end module
