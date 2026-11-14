module exclusion_reset_mod
  use iso_c_binding,only:c_ptr,c_null_ptr,c_associated,c_f_pointer,c_loc,c_int,c_size_t,c_int8_t
  implicit none
  private
  public::exclusion_reset_enabled,exclusion_reset_verify,exclusion_restore_masks
  type(c_ptr),save::previous_one=c_null_ptr,previous_inter=c_null_ptr
  integer,save::previous_ld=0,previous_local=0,previous_near=0
  interface
    function host_memset(address,value,bytes) bind(C,name='memset') result(destination)
      import c_ptr,c_int,c_size_t
      type(c_ptr),value::address
      integer(c_int),value::value
      integer(c_size_t),value::bytes
      type(c_ptr)::destination
    end function
  end interface
contains
  logical function flag(name)
    character(*),intent(in)::name
    character(16)::value
    integer::status
    call get_environment_variable(name,value,status=status)
    flag=.false.
    if(status==1)return
    if(status/=0)error stop 'Truncated exclusion reset flag'
    if(trim(value)/='0'.and.trim(value)/='1')error stop 'Invalid exclusion reset flag'
    flag=trim(value)=='1'
  end function
  logical function exclusion_reset_enabled()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_EXCLUSION_RESET');initialized=.true.
    end if
    exclusion_reset_enabled=selected
  end function
  logical function exclusion_reset_verify()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_EXCLUSION_VERIFY');initialized=.true.
    end if
    exclusion_reset_verify=selected
  end function
  subroutine exclusion_restore_masks(ld,local,near,one_address,inter_address,records,used)
    integer,intent(in)::ld,local,near
    type(c_ptr),value::one_address,inter_address
    integer,optional,intent(in)::records(:,:,:),used(:,:)
    integer(c_int8_t),pointer,contiguous::one(:,:,:),inter(:,:,:)
    integer::p,row,thread,entry,column,plane,restored
    logical::full
    type(c_ptr)::destination
    call c_f_pointer(one_address,one,[ld,ld,local])
    call c_f_pointer(inter_address,inter,[ld,ld,near])
    full=ld/=previous_ld.or.local/=previous_local.or.near/=previous_near.or. &
      .not.c_associated(previous_one,one_address).or..not.c_associated(previous_inter,inter_address)
    restored=0
    if(full)then
      !$omp parallel do private(row,destination) schedule(static)
      do p=1,local+near
        if(p<=local)then
          destination=host_memset(c_loc(one(1,1,p)),1_c_int,int(ld,c_size_t)*ld)
          do row=1,ld
            one(row,row,p)=0_c_int8_t
          end do
        else
          destination=host_memset(c_loc(inter(1,1,p-local)),1_c_int,int(ld,c_size_t)*ld)
        end if
      end do
      !$omp end parallel do
    else
      if(.not.present(records).or..not.present(used))error stop 'Missing previous exclusion records'
      if(size(records,3)/=size(used,2))error stop 'Exclusion thread-history mismatch'
      restored=sum(used(1,:))
      !$omp parallel do private(entry,row,column,plane) schedule(static)
      do thread=1,size(used,2)
        if(used(1,thread)<0.or.used(1,thread)>size(records,2))error stop 'Invalid exclusion history size'
        do entry=1,used(1,thread)
          row=records(1,entry,thread);column=records(2,entry,thread);plane=records(3,entry,thread)
          if(row<1.or.row>ld.or.column<1.or.column>ld)error stop 'Invalid old exclusion index'
          if(plane<0)then
            plane=-plane
            if(plane>local.or.row==column)error stop 'Invalid old intra exclusion'
            ! Duplicate records may cross threads; atomic byte stores avoid a data race.
            !$omp atomic write
            one(row,column,plane)=1_c_int8_t
          else
            if(plane<1.or.plane>near)error stop 'Invalid old inter exclusion'
            !$omp atomic write
            inter(row,column,plane)=1_c_int8_t
          end if
        end do
      end do
      !$omp end parallel do
    end if
    previous_ld=ld;previous_local=local;previous_near=near
    previous_one=one_address;previous_inter=inter_address
    if(flag('GENESIS_EXCLUSION_STATS')) &
      write(*,'(a,l1,4i12)') 'GENESIS_EXCLUSION_RESET: full,ld,local,near,records ',full,ld,local,near,restored
  end subroutine
end module
