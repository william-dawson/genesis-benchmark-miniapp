module exclusion_oracle_mod
  use constants_mod
  use iso_c_binding,only:c_ptr,c_loc,c_int,c_size_t
  use sp_enefunc_str_mod
  implicit none
  private
  public::exclusion_save,exclusion_check
  integer,allocatable,target,save::counts(:),counts14(:),lists(:,:,:),lists14(:,:,:),scales(:,:)
  real(wp),allocatable,target,save::qq(:,:),lj(:,:)
  integer(1),allocatable,target,save::one(:,:,:),inter(:,:,:)
  interface
    function copy(a,b,n) bind(C,name='memcpy') result(destination)
      import c_ptr,c_size_t
      type(c_ptr),value::a,b
      integer(c_size_t),value::n
      type(c_ptr)::destination
    end function
    integer(c_int) function compare(a,b,n) bind(C,name='memcmp')
      import c_ptr,c_int,c_size_t
      type(c_ptr),value::a,b
      integer(c_size_t),value::n
    end function
  end interface
contains
  subroutine exclusion_save(enefunc)
    type(s_enefunc),target,intent(in)::enefunc
    integer::p,local,near,ld
    type(c_ptr)::destination
    counts=enefunc%num_nonb_excl;counts14=enefunc%num_nb14_calc
    lists=enefunc%nonb_excl_list;lists14=enefunc%nb14_calc_list;scales=enefunc%sc_calc_list
    qq=enefunc%nb14_qq_scale;lj=enefunc%nb14_lj_scale
    ld=size(enefunc%exclusion_mask1,1);local=size(enefunc%exclusion_mask1,3);near=size(enefunc%exclusion_mask,3)
    if(allocated(one))then
      if(any(shape(one)/=shape(enefunc%exclusion_mask1)))deallocate(one)
    end if
    if(allocated(inter))then
      if(any(shape(inter)/=shape(enefunc%exclusion_mask)))deallocate(inter)
    end if
    if(.not.allocated(one))allocate(one(ld,ld,local))
    if(.not.allocated(inter))allocate(inter(ld,ld,near))
    !$omp parallel do private(destination) schedule(static)
    do p=1,local+near
      if(p<=local)then
        destination=copy(c_loc(one(1,1,p)),c_loc(enefunc%exclusion_mask1(1,1,p)),int(ld,c_size_t)*ld)
      else
        destination=copy(c_loc(inter(1,1,p-local)),c_loc(enefunc%exclusion_mask(1,1,p-local)),int(ld,c_size_t)*ld)
      end if
    end do
    !$omp end parallel do
  end subroutine
  subroutine equal(a,b,n)
    type(c_ptr),value::a,b
    integer(c_size_t),intent(in)::n
    if(compare(a,b,n)/=0)error stop 'Exclusion rebuild differs from full reference'
  end subroutine
  subroutine exclusion_check(enefunc,ld,local,near,natom,pairs)
    type(s_enefunc),target,intent(in)::enefunc
    integer,intent(in)::ld,local,near,natom(*)
    integer(int2),intent(in)::pairs(2,*)
    integer::cell,plane,n,k
    integer(c_size_t)::bytes
    if(any(counts/=enefunc%num_nonb_excl).or.any(counts14/=enefunc%num_nb14_calc)) &
      error stop 'Exclusion/1-4 count mismatch'
    !$omp parallel do private(n,k,bytes) schedule(static)
    do cell=1,local
      n=counts(cell)
      if(n>0)call equal(c_loc(lists(1,1,cell)),c_loc(enefunc%nonb_excl_list(1,1,cell)),int(16*n,c_size_t))
      n=counts14(cell)
      if(n>0)then
        call equal(c_loc(lists14(1,1,cell)),c_loc(enefunc%nb14_calc_list(1,1,cell)),int(16*n,c_size_t))
        call equal(c_loc(scales(1,cell)),c_loc(enefunc%sc_calc_list(1,cell)),int(4*n,c_size_t))
        bytes=int(storage_size(qq)/8*n,c_size_t)
        call equal(c_loc(qq(1,cell)),c_loc(enefunc%nb14_qq_scale(1,cell)),bytes)
        call equal(c_loc(lj(1,cell)),c_loc(enefunc%nb14_lj_scale(1,cell)),bytes)
      end if
      k=natom(cell)
      if(k>0)call equal(c_loc(one(1,1,cell)),c_loc(enefunc%exclusion_mask1(1,1,cell)), &
                        int(k-1,c_size_t)*ld+k)
    end do
    !$omp end parallel do
    !$omp parallel do private(k) schedule(static)
    do plane=1,near
      k=max(natom(pairs(1,plane)),natom(pairs(2,plane)))
      if(k>0)call equal(c_loc(inter(1,1,plane)),c_loc(enefunc%exclusion_mask(1,1,plane)), &
                        int(k-1,c_size_t)*ld+k)
    end do
    !$omp end parallel do
    write(*,'(a,2i12)') 'GENESIS_EXCLUSION_VERIFY: exact masks and ordered lists pass ',sum(counts),sum(counts14)
  end subroutine
end module
