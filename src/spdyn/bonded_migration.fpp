module bonded_migration_mod
  use constants_mod
  use iso_c_binding,only:c_ptr,c_int,c_int16_t,c_double,c_size_t,c_loc,c_f_pointer
  use sp_domain_str_mod
  use sp_enefunc_str_mod
  use mpi_parallel_mod
  implicit none
  private
  public::bonded_outgoing,bonded_check,bonded_verify,bonded_update_ids,bonded_check_ids
  integer,allocatable,save::starts(:),saved_add(:),saved_exit(:),add_start(:),exit_start(:)
  integer,allocatable,target,save::saved_indices(:),saved_int(:,:)
  real(wp),allocatable,target,save::saved_real(:,:)
  type topology_state
    integer,allocatable::ids(:,:),incoming(:,:)
    integer,allocatable::initial_int(:,:)
    real(wp),allocatable::initial_real(:,:)
    integer::total=0
  end type
  type(topology_state),target,save::topology(3)
  integer,pointer,contiguous,save::batch_exits(:),batch_ids(:),batch_int(:,:),source_ranges(:,:),destination_ranges(:,:)
  real(wp),pointer,contiguous,save::batch_real(:,:)
  interface
    function staging(n) bind(C,name='genesis_bonded_staging') result(p)
      import c_ptr,c_int
      integer(c_int),value::n
      type(c_ptr)::p
    end function
    subroutine cache_c(kind,n,integers,reals) bind(C,name='genesis_bonded_cache')
      import c_int,c_double
      integer(c_int),value::kind,n
      integer(c_int),intent(in)::integers(*)
      real(c_double),intent(in)::reals(*)
    end subroutine
    subroutine batch_c(n,map,owner,atoms,cells,pitch,moved,exits,ids,integers,reals,src,dst) &
        bind(C,name='genesis_bonded_batch')
      import c_ptr,c_int,c_int16_t
      integer(c_int),value::n,atoms,cells,pitch
      integer(c_int16_t),intent(in)::map(*),owner(*)
      integer(c_int),intent(out)::moved
      type(c_ptr),intent(out)::exits,ids,integers,reals,src,dst
    end subroutine
    integer(c_int) function compare(a,b,n) bind(C,name='memcmp')
      import c_ptr,c_int,c_size_t
      type(c_ptr),value::a,b
      integer(c_size_t),value::n
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
    if(status/=0)error stop 'Truncated bonded migration flag'
    if(trim(value)/='0'.and.trim(value)/='1')error stop 'Invalid bonded migration flag'
    flag=trim(value)=='1'
  end function
  logical function enabled()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_BONDED_MIGRATION_GPU');initialized=.true.
    end if
    enabled=selected
  end function
  logical function bonded_verify()
    logical,save::initialized=.false.,selected=.false.
    if(.not.initialized)then
      selected=flag('GENESIS_BONDED_MIGRATION_VERIFY');initialized=.true.
    end if
    bonded_verify=selected
  end function
  real(8) function now()
    integer(8)::ticks,rate
    call system_clock(ticks,rate)
    now=real(ticks,8)/real(rate,8)
  end function
  subroutine fields(kind,e,counts,adds,exits,indices,integer_buffer,real_buffer,lists)
    integer,intent(in)::kind
    type(s_enefunc),target,intent(inout)::e
    integer,pointer,intent(out)::counts(:),adds(:),exits(:),indices(:,:),integer_buffer(:,:,:),lists(:,:,:)
    real(wp),pointer,intent(out)::real_buffer(:,:,:)
    select case(kind)
    case(1)
      counts=>e%num_bond;adds=>e%bond_add;exits=>e%bond_exit;indices=>e%bond_exit_index
      integer_buffer=>e%buf_bond_integer;real_buffer=>e%buf_bond_real;lists=>e%bond_list
    case(2)
      counts=>e%num_angle;adds=>e%angle_add;exits=>e%angle_exit;indices=>e%angle_exit_index
      integer_buffer=>e%buf_angle_integer;real_buffer=>e%buf_angle_real;lists=>e%angle_list
    case(3)
      counts=>e%num_dihedral;adds=>e%dihed_add;exits=>e%dihed_exit;indices=>e%dihed_exit_index
      integer_buffer=>e%buf_dihed_integer;real_buffer=>e%buf_dihed_real;lists=>e%dihe_list
    case default
      error stop 'Invalid bonded migration kind'
    end select
  end subroutine
  subroutine payload(kind,e,cell,slot,iv,rv)
    integer,intent(in)::kind,cell,slot
    type(s_enefunc),intent(in)::e
    integer,intent(out)::iv(9)
    real(wp),intent(out)::rv(4)
    iv=0;rv=0
    select case(kind)
    case(1)
      iv(1:4)=[e%bond_list(:,slot,cell),int(e%bond_kind(slot,cell)),e%bond_pbc(slot,cell)]
      rv(1:2)=[e%bond_dist_min(slot,cell),e%bond_force_const(slot,cell)]
    case(2)
      iv(1:7)=[e%angle_list(:,slot,cell),int(e%angle_kind(slot,cell)),e%angle_pbc(:,slot,cell)]
      rv=[e%angle_theta_min(slot,cell),e%angle_force_const(slot,cell), &
          e%urey_rmin(slot,cell),e%urey_force_const(slot,cell)]
    case(3)
      iv=[e%dihe_list(:,slot,cell),e%dihe_periodicity(slot,cell),int(e%dihe_kind(slot,cell)),e%dihe_pbc(:,slot,cell)]
      rv(1:2)=[e%dihe_phase(slot,cell),e%dihe_force_const(slot,cell)]
    end select
  end subroutine
  subroutine prepare(domain,e,prep_us,driver_us)
    type(s_domain),target,intent(in)::domain
    type(s_enefunc),target,intent(inout)::e
    real(8),intent(out)::prep_us,driver_us
    integer,pointer::counts(:),adds(:),exits(:),indices(:,:),integer_buffer(:,:,:),lists(:,:,:)
    real(wp),pointer::real_buffer(:,:,:)
    integer,pointer,contiguous::records(:,:)
    integer::kind,cell,slot,index,total,cells,n,moved,offset
    type(c_ptr)::exit_p,id_p,int_p,real_p,src_p,dst_p
    real(8)::start,ready_time
    start=now();cells=domain%num_cell_local
    if(allocated(starts))then
      if(size(starts)/=cells)error stop 'Bonded cell layout changed'
    else
      allocate(starts(cells))
    end if
    total=sum(e%num_bond(1:cells))+sum(e%num_angle(1:cells))+sum(e%num_dihedral(1:cells))
    call c_f_pointer(staging(max(1,total)),records,[4,max(1,total)])
    offset=0
    do kind=1,3
      call fields(kind,e,counts,adds,exits,indices,integer_buffer,real_buffer,lists)
      n=0
      do cell=1,cells
        starts(cell)=n;n=n+counts(cell)
      end do
      if(.not.allocated(topology(kind)%ids))then
        topology(kind)%total=n
        allocate(topology(kind)%ids(size(lists,2),cells),topology(kind)%incoming(size(integer_buffer,2),cells), &
                 topology(kind)%initial_int(9,n),topology(kind)%initial_real(4,n))
        topology(kind)%ids=0;topology(kind)%incoming=0
        !$omp parallel do private(slot,index) schedule(static)
        do cell=1,cells
          do slot=1,counts(cell)
            index=starts(cell)+slot;topology(kind)%ids(slot,cell)=index
            call payload(kind,e,cell,slot,topology(kind)%initial_int(:,index),topology(kind)%initial_real(:,index))
          end do
        end do
        !$omp end parallel do
        call cache_c(kind,n,topology(kind)%initial_int,topology(kind)%initial_real)
      end if
      if(n/=topology(kind)%total.or.size(lists,2)/=size(topology(kind)%ids,1)) &
        error stop 'Bonded topology changed'
      !$omp parallel do private(slot,index) schedule(static)
      do cell=1,cells
        do slot=1,counts(cell)
          index=offset+starts(cell)+slot
          records(1,index)=topology(kind)%ids(slot,cell);records(2,index)=0
          records(3,index)=cell+ishft(kind-1,16);records(4,index)=slot
        end do
      end do
      !$omp end parallel do
      offset=offset+n
    end do
    ready_time=now();prep_us=(ready_time-start)*1d6
    call batch_c(total,domain%id_g2l,domain%cell_pair,size(domain%id_g2l,2),cells,size(domain%cell_pair,1), &
                 moved,exit_p,id_p,int_p,real_p,src_p,dst_p)
    call c_f_pointer(exit_p,batch_exits,[max(1,moved)])
    call c_f_pointer(id_p,batch_ids,[max(1,moved)])
    call c_f_pointer(int_p,batch_int,[9,max(1,moved)])
    call c_f_pointer(real_p,batch_real,[4,max(1,moved)])
    call c_f_pointer(src_p,source_ranges,[2,3*cells])
    call c_f_pointer(dst_p,destination_ranges,[2,3*cells])
    driver_us=(now()-ready_time)*1d6
  end subroutine
  logical function bonded_outgoing(kind,domain,e)
    integer,intent(in)::kind
    type(s_domain),target,intent(in)::domain
    type(s_enefunc),target,intent(inout)::e
    integer,pointer::counts(:),adds(:),exits(:),indices(:,:),integer_buffer(:,:,:),lists(:,:,:)
    real(wp),pointer::real_buffer(:,:,:)
    integer::cell,cells,b,first,last,n,ni,nr
    real(8)::prep_us,driver_us,start
    logical,save::timing_initialized=.false.,timing=.false.
    bonded_outgoing=.false.
    if(.not.enabled())return
    if(nproc_country/=1.or.domain%num_cell_boundary/=0.or.e%forcefield/=ForcefieldCHARMM.or. &
       domain%nonbond_kernel/=NBK_GPU.or.storage_size(1.0_wp)/=64) &
      error stop 'GPU bonded migration requires one-rank CHARMM double GPU and no halos'
    cells=domain%num_cell_local
    if(cells>=32768.or.size(domain%cell_pair,1)<cells.or.size(domain%cell_pair,1)/=size(domain%cell_pair,2)) &
      error stop 'Unexpected bonded cell-owner layout'
    if(.not.timing_initialized)then
      timing=flag('GENESIS_BONDED_MIGRATION_TIMING');timing_initialized=.true.
    end if
    prep_us=0;driver_us=0
    if(kind==1)call prepare(domain,e,prep_us,driver_us)
    call fields(kind,e,counts,adds,exits,indices,integer_buffer,real_buffer,lists)
    ni=size(integer_buffer,1);nr=size(real_buffer,1);start=now()
    !$omp parallel do private(b,first,last,n) schedule(static)
    do cell=1,cells
      b=cell+(kind-1)*cells
      first=source_ranges(1,b);last=source_ranges(2,b);n=0
      if(first>0)n=last-first+1
      if(n>size(indices,1))error stop 'Bonded exit capacity exceeded'
      exits(cell)=n
      if(n>0)indices(1:n,cell)=batch_exits(first:last)
      first=destination_ranges(1,b);last=destination_ranges(2,b);n=0
      if(first>0)n=last-first+1
      if(n>size(integer_buffer,2))error stop 'Bonded incoming capacity exceeded'
      adds(cell)=n
      if(n>0)then
        integer_buffer(:,1:n,cell)=batch_int(1:ni,first:last)
        real_buffer(:,1:n,cell)=batch_real(1:nr,first:last)
        topology(kind)%incoming(1:n,cell)=batch_ids(first:last)
      end if
    end do
    !$omp end parallel do
    if(timing)write(*,'(a,3i12,3es20.10)') 'GENESIS_BONDED_MIGRATION_HOST: ', &
      kind,topology(kind)%total,sum(adds(1:cells)),prep_us,driver_us,(now()-start)*1d6
    if(bonded_verify())call save(cells,adds,exits,indices,integer_buffer,real_buffer)
    bonded_outgoing=.true.
  end function
  subroutine bonded_update_ids(kind,e)
    integer,intent(in)::kind
    type(s_enefunc),target,intent(inout)::e
    integer,pointer::counts(:),adds(:),exits(:),indices(:,:),integers(:,:,:),lists(:,:,:)
    real(wp),pointer::reals(:,:,:)
    integer::cell,k,j,ix,kx,n
    logical::insert
    if(.not.enabled())return
    call fields(kind,e,counts,adds,exits,indices,integers,reals,lists)
    if(.not.allocated(topology(kind)%ids))error stop 'Missing bonded topology IDs'
    !$omp parallel do private(k,j,ix,kx,n,insert) schedule(static)
    do cell=1,size(topology(kind)%ids,2)
      n=counts(cell)
      if(adds(cell)>=exits(cell))then
        do k=1,exits(cell)
          topology(kind)%ids(indices(k,cell),cell)=topology(kind)%incoming(k,cell)
        end do
        do k=exits(cell)+1,adds(cell)
          topology(kind)%ids(k+n-exits(cell),cell)=topology(kind)%incoming(k,cell)
        end do
      else
        do k=1,adds(cell)
          topology(kind)%ids(indices(k,cell),cell)=topology(kind)%incoming(k,cell)
        end do
        j=0;ix=n;k=adds(cell)+1
        do while(j<exits(cell)-adds(cell))
          insert=.true.
          do kx=k,exits(cell)
            if(ix==indices(kx,cell))then
              insert=.false.;j=j+1;exit
            end if
          end do
          if(insert)then
            topology(kind)%ids(indices(k,cell),cell)=topology(kind)%ids(ix,cell)
            j=j+1;k=k+1
          end if
          ix=ix-1
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine
  subroutine bonded_check_ids(kind,e)
    integer,intent(in)::kind
    type(s_enefunc),target,intent(inout)::e
    integer,pointer::counts(:),adds(:),exits(:),indices(:,:),integers(:,:,:),lists(:,:,:)
    real(wp),pointer::reals(:,:,:)
    integer,allocatable::seen(:)
    integer,target::iv(9)
    real(wp),target::rv(4)
    integer::cell,slot,id,ni,nr
    if(.not.enabled().or..not.bonded_verify())return
    call fields(kind,e,counts,adds,exits,indices,integers,reals,lists)
    ni=size(integers,1);nr=size(reals,1)
    allocate(seen(topology(kind)%total));seen=0
    do cell=1,size(topology(kind)%ids,2)
      do slot=1,counts(cell)
        id=topology(kind)%ids(slot,cell)
        if(id<1.or.id>size(seen))error stop 'Invalid bonded topology ID'
        seen(id)=seen(id)+1
        call payload(kind,e,cell,slot,iv,rv)
        if(any(iv(1:ni)/=topology(kind)%initial_int(1:ni,id)))error stop 'Bonded topology integer mismatch'
        if(compare(c_loc(rv(1)),c_loc(topology(kind)%initial_real(1,id)),int(8*nr,c_size_t))/=0) &
          error stop 'Bonded topology real mismatch'
      end do
    end do
    if(any(seen/=1))error stop 'Bonded topology ID conservation mismatch'
    write(*,'(a,2i12)') 'GENESIS_BONDED_MIGRATION_VERIFY: exact incoming topology IDs pass ',kind,size(seen)
  end subroutine
  subroutine save(cells,adds,exits,indices,integers,reals)
    integer,intent(in)::cells,adds(:),exits(:),indices(:,:),integers(:,:,:)
    real(wp),intent(in)::reals(:,:,:)
    integer::cell,a,b
    if(allocated(saved_add))deallocate(saved_add,saved_exit,add_start,exit_start,saved_indices,saved_int,saved_real)
    allocate(saved_add(cells),saved_exit(cells),add_start(cells),exit_start(cells))
    saved_add=adds(1:cells);saved_exit=exits(1:cells)
    allocate(saved_indices(sum(saved_exit)),saved_int(size(integers,1),sum(saved_add)), &
             saved_real(size(reals,1),sum(saved_add)))
    a=0;b=0
    do cell=1,cells
      add_start(cell)=a;exit_start(cell)=b
      saved_indices(b+1:b+exits(cell))=indices(1:exits(cell),cell)
      saved_int(:,a+1:a+adds(cell))=integers(:,1:adds(cell),cell)
      saved_real(:,a+1:a+adds(cell))=reals(:,1:adds(cell),cell)
      a=a+adds(cell);b=b+exits(cell)
    end do
  end subroutine
  subroutine bonded_check(kind,e)
    integer,intent(in)::kind
    type(s_enefunc),target,intent(inout)::e
    integer,pointer::counts(:),adds(:),exits(:),indices(:,:),integers(:,:,:),lists(:,:,:)
    real(wp),pointer::reals(:,:,:)
    integer::cell,cells,a,b,n
    if(.not.enabled().or..not.bonded_verify())return
    call fields(kind,e,counts,adds,exits,indices,integers,reals,lists)
    cells=size(saved_add)
    if(any(saved_add/=adds(1:cells)).or.any(saved_exit/=exits(1:cells)))error stop 'Bonded migration count mismatch'
    do cell=1,cells
      a=add_start(cell);b=exit_start(cell);n=exits(cell)
      if(n>0)then
        if(compare(c_loc(saved_indices(b+1)),c_loc(indices(1,cell)),int(4*n,c_size_t))/=0) &
          error stop 'Bonded migration exit order mismatch'
      end if
      n=adds(cell)
      if(n>0)then
        if(compare(c_loc(saved_int(1,a+1)),c_loc(integers(1,1,cell)),int(4*size(integers,1)*n,c_size_t))/=0) &
          error stop 'Bonded migration integer payload mismatch'
        if(compare(c_loc(saved_real(1,a+1)),c_loc(reals(1,1,cell)),int(storage_size(reals)/8*size(reals,1)*n,c_size_t))/=0) &
          error stop 'Bonded migration real payload mismatch'
      end if
    end do
    write(*,'(a,3i12)') 'GENESIS_BONDED_MIGRATION_VERIFY: exact counts exits and payload pass ',kind, &
      sum(saved_exit),sum(saved_add)
  end subroutine
end module
