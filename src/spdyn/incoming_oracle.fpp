module incoming_oracle_mod
  use constants_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use iso_c_binding,only:c_ptr,c_loc,c_size_t,c_int
  implicit none
  private
  public::incoming_mode,incoming_verify,incoming_save,incoming_check
  real(wip),allocatable,target,save::x(:,:,:),v(:,:,:),m(:,:),dist(:,:,:,:)
  real(wp),allocatable,target,save::q(:,:)
  integer,allocatable,target,save::ids(:,:),sol(:,:),cls(:,:),water(:,:,:),groups(:,:,:,:)
  integer,allocatable,save::nw(:),na(:),ns(:),ng(:,:)
  integer(int2),allocatable,target,save::reverse(:,:)
  interface
    integer(c_int) function compare(a,b,n) bind(C,name='memcmp')
      import c_ptr,c_int,c_size_t
      type(c_ptr),value::a,b
      integer(c_size_t),value::n
    end function
  end interface
contains
  integer function incoming_mode()
    character(16)::value
    integer::status,code
    integer,save::mode=-1
    if(mode<0)then
      mode=0
      call get_environment_variable('GENESIS_MIGRATION_HOST',value,status=status)
      if(status==0)then
        read(value,*,iostat=code)mode
        if(code/=0.or.mode<0.or.mode>3)error stop 'Invalid GENESIS_MIGRATION_HOST'
      end if
    end if
    incoming_mode=mode
  end function
  logical function incoming_verify()
    character(16)::value
    integer::status
    call get_environment_variable('GENESIS_MIGRATION_VERIFY',value,status=status)
    incoming_verify=status==0.and.trim(value)=='1'
  end function
  subroutine incoming_save(domain,constraints,kind)
    type(s_domain),intent(in)::domain
    type(s_constraints),optional,intent(in)::constraints
    integer,intent(in)::kind
    x=domain%coord;v=domain%velocity;m=domain%mass;q=domain%charge
    ids=domain%id_l2g;cls=domain%atom_cls_no
    if(kind==1)then
      nw=domain%num_water;na=domain%num_atom;water=domain%water_list
    else
      ns=domain%num_solute;sol=domain%id_l2g_solute;reverse=domain%id_g2l
      ng=constraints%HGr_local;groups=constraints%HGr_bond_list;dist=constraints%HGr_bond_dist
    end if
  end subroutine
  subroutine equal(a,b,n)
    type(c_ptr),value::a,b
    integer(c_size_t),intent(in)::n
    if(compare(a,b,n)/=0)error stop 'Incoming migration differs from serial reference'
  end subroutine
  subroutine incoming_check(domain,constraints,kind)
    type(s_domain),target,intent(in)::domain
    type(s_constraints),optional,target,intent(in)::constraints
    integer,intent(in)::kind
    integer::cell,first,last,n,j,k,id,sites
    integer(c_size_t)::bytes
    do cell=1,domain%num_cell_local
      if(kind==1)then
        if(nw(cell)/=domain%num_water(cell).or.na(cell)/=domain%num_atom(cell))error stop 'Water count mismatch'
        first=domain%num_solute(cell)+1;last=domain%num_atom(cell)
        if(nw(cell)>0)call equal(c_loc(water(1,1,cell)),c_loc(domain%water_list(1,1,cell)),int(12*nw(cell),c_size_t))
      else
        if(ns(cell)/=domain%num_solute(cell))error stop 'Solute count mismatch'
        first=constraints%No_HGr(cell)+1;last=domain%num_solute(cell)
        do j=1,constraints%connect
          n=constraints%HGr_local(j,cell)
          if(ng(j,cell)/=n)error stop 'Constrained-group count mismatch'
          do k=1,n
            sites=j+1
            call equal(c_loc(groups(1,k,j,cell)),c_loc(constraints%HGr_bond_list(1,k,j,cell)),int(4*sites,c_size_t))
            call equal(c_loc(dist(1,k,j,cell)),c_loc(constraints%HGr_bond_dist(1,k,j,cell)), &
              int(storage_size(dist)/8*sites,c_size_t))
          end do
        end do
      end if
      n=last-first+1
      if(n<=0)cycle
      bytes=int(storage_size(x)/8*n,c_size_t)
      call equal(c_loc(x(1,first,cell)),c_loc(domain%coord(1,first,cell)),3*bytes)
      call equal(c_loc(v(1,first,cell)),c_loc(domain%velocity(1,first,cell)),3*bytes)
      call equal(c_loc(m(first,cell)),c_loc(domain%mass(first,cell)),bytes)
      call equal(c_loc(q(first,cell)),c_loc(domain%charge(first,cell)),int(storage_size(q)/8*n,c_size_t))
      call equal(c_loc(ids(first,cell)),c_loc(domain%id_l2g(first,cell)),int(4*n,c_size_t))
      call equal(c_loc(cls(first,cell)),c_loc(domain%atom_cls_no(first,cell)),int(4*n,c_size_t))
      if(kind==2)then
        call equal(c_loc(sol(first,cell)),c_loc(domain%id_l2g_solute(first,cell)),int(4*n,c_size_t))
        do k=first,last
          id=domain%id_l2g_solute(k,cell)
          call equal(c_loc(reverse(1,id)),c_loc(domain%id_g2l(1,id)),4_c_size_t)
        end do
      end if
    end do
    write(*,'(a,i0)') 'GENESIS_MIGRATION_VERIFY: exact incoming pass ',kind
  end subroutine
end module
