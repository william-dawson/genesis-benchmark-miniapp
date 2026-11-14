#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module genesis_pair_sparse_mod
  use exclusion_reset_mod,only:exclusion_restore_masks
  use iso_c_binding, only: c_ptr, c_loc, c_int, c_f_pointer
#ifdef OMP
  use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
  implicit none
  private
  public :: pair_sparse_enabled, pair_sparse_reset, pair_sparse_record, pair_sparse_submit
  public :: pair_sparse_restore
  integer, parameter :: capacity=65536
  integer, allocatable, target, save :: records(:,:,:), compact(:,:), inverse(:)
  integer, allocatable, save :: used(:,:)
  interface
    function gpu_sparse_staging(count) bind(c,name='gpu_sparse_staging_') result(address)
      import :: c_ptr, c_int
      integer(c_int), value :: count
      type(c_ptr) :: address
    end function
    subroutine gpu_sparse_mask(entries, mapping, count, width, pitch, local, near) bind(c,name='gpu_sparse_mask_')
      import :: c_ptr, c_int
      type(c_ptr), value :: entries, mapping
      integer(c_int), value :: count, width, pitch, local, near
    end subroutine
  end interface
contains
  subroutine pair_sparse_restore(ld,local,near,one,inter)
    integer,intent(in)::ld,local,near
    type(c_ptr),value::one,inter
    if(allocated(records))then
      call exclusion_restore_masks(ld,local,near,one,inter,records,used)
    else
      call exclusion_restore_masks(ld,local,near,one,inter)
    end if
  end subroutine

  logical function pair_sparse_enabled()
    integer, save :: selected=-1
    integer :: status, code
    character(16) :: value
    if (selected<0) then
      selected=0
      call get_environment_variable('GENESIS_PAIR_GPU',value,status=status)
      if (status==0) then
        read(value,*,iostat=code) selected
        if (code/=0) error stop 'Invalid GENESIS_PAIR_GPU'
        if (selected<0 .or. selected>1) error stop 'GENESIS_PAIR_GPU must be 0 or 1'
      end if
    end if
    pair_sparse_enabled=selected==1
  end function

  subroutine pair_sparse_reset()
    integer :: threads, stride, status, code, padding
    character(16) :: value
    threads=1
#ifdef OMP
    threads=omp_get_max_threads()
#endif
    padding=0
    call get_environment_variable('GENESIS_PAIR_RECORD_PAD',value,status=status)
    if (status==0) then
      read(value,*,iostat=code) padding
      if (code/=0) error stop 'Invalid GENESIS_PAIR_RECORD_PAD'
      if (padding<0 .or. padding>1) error stop 'GENESIS_PAIR_RECORD_PAD must be 0 or 1'
    end if
    stride=1
    if (padding==1) stride=32
    if (allocated(used)) then
      if (size(used,2)/=threads .or. size(used,1)/=stride) deallocate(records,used)
    end if
    if (.not. allocated(used)) allocate(records(4,capacity,threads),used(stride,threads))
    used=0
  end subroutine

  subroutine pair_sparse_record(row,column,plane,intra)
    integer, intent(in) :: row,column,plane
    logical, intent(in) :: intra
    integer :: thread, entry
    if (intra .and. row==column) return
    thread=1
#ifdef OMP
    thread=omp_get_thread_num()+1
#endif
    entry=used(1,thread)+1
    if (entry>capacity) error stop 'Sparse mask record capacity exceeded'
    used(1,thread)=entry
    records(1,entry,thread)=row
    records(2,entry,thread)=column
    records(3,entry,thread)=plane
    if (intra) records(3,entry,thread)=-plane
    records(4,entry,thread)=0
  end subroutine

  subroutine pair_sparse_submit(width,pitch,local,near,pairs)
    integer, intent(in) :: width,pitch,local,near,pairs(2,*)
    integer :: total,thread,offset,plane,status,code,pinned
    integer, allocatable :: starts(:)
    integer, pointer, contiguous :: target_records(:,:)
    type(c_ptr) :: address
    character(16) :: value
    if (.not. allocated(used)) error stop 'Missing sparse exclusion records'
    total=sum(used)
    pinned=0
    call get_environment_variable('GENESIS_PAIR_PINNED',value,status=status)
    if (status==0) then
      read(value,*,iostat=code) pinned
      if (code/=0) error stop 'Invalid GENESIS_PAIR_PINNED'
      if (pinned<0 .or. pinned>1) error stop 'GENESIS_PAIR_PINNED must be 0 or 1'
    end if
    if (pinned==1) then
#ifdef USE_GPU
      address=gpu_sparse_staging(total)
      call c_f_pointer(address,target_records,[4,max(1,total)])
#else
      error stop 'Pinned sparse mask requires a CUDA build'
#endif
    else
      if (allocated(compact)) then
        if (size(compact,2)<max(1,total)) deallocate(compact)
      end if
      if (.not. allocated(compact)) allocate(compact(4,max(1,total+total/10)))
      target_records=>compact
    end if
    allocate(starts(size(used,2)))
    offset=0
    do thread=1,size(used,2)
      starts(thread)=offset
      offset=offset+used(1,thread)
    end do
    !$omp parallel do schedule(static) private(offset)
    do thread=1,size(used,2)
      offset=starts(thread)
      target_records(:,offset+1:offset+used(1,thread))=records(:,1:used(1,thread),thread)
    end do
    !$omp end parallel do
    if (allocated(inverse)) then
      if (size(inverse)/=local) deallocate(inverse)
    end if
    if (.not. allocated(inverse)) allocate(inverse(local))
    inverse=0
    do plane=1,local
      if (pairs(1,plane)<1 .or. pairs(1,plane)>local) error stop 'Invalid local mask cell'
      inverse(pairs(1,plane))=plane
    end do
    if (any(inverse==0)) error stop 'Missing local mask plane'
#ifdef USE_GPU
    call gpu_sparse_mask(c_loc(target_records),c_loc(inverse),total,width,pitch,local,near)
#else
    error stop 'Sparse GPU mask requires a CUDA build'
#endif
  end subroutine
end module
