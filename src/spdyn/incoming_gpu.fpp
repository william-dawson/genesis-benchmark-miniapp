#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module incoming_gpu_mod
  use iso_c_binding
  use constants_mod
  use mpi_parallel_mod
  use sp_domain_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod,only:TIP1,TIP2,TIP3,TIP4
  use sp_constraints_str_mod
  use incoming_oracle_mod,only:incoming_verify
  implicit none
  private
  public::incoming_gpu_enabled,incoming_gpu_verify,incoming_atom_gpu,incoming_water_gpu,incoming_group_gpu
  public::outgoing_atom_gpu,outgoing_atom_verify,outgoing_atom_check
  public::incoming_atom_backup,incoming_atom_restore
  public::incoming_boundary_gpu,incoming_boundary_check,incoming_boundary_suspend
  type,bind(C)::particle_view
    type(c_ptr)::x,v,q,m,classes,ids,solute,reverse,ns,nw,na,noh,groups,distance,hlocal,water
    integer(c_int)::width,cells,water_height,water_capacity,connect,group_height,group_capacity
    integer(c_int)::reverse_capacity,dist_height,dist_capacity
  end type
  type,bind(C)::records
    type(c_ptr)::move_real,stay_real,move_int,stay_int,move_count,stay_count
    integer(c_int)::move_real_rows,move_real_capacity,stay_real_rows,stay_real_capacity
    integer(c_int)::move_int_rows,move_int_capacity,stay_int_rows,stay_int_capacity
  end type
  type,bind(C)::atom_records
    type(c_ptr)::reals,integers,adds,exits,indices
    integer(c_int)::real_rows,real_capacity,int_rows,int_capacity,index_capacity
  end type
  real(wip),allocatable,save::old_x(:,:,:),old_v(:,:,:),old_m(:,:)
  real(wp),allocatable,save::old_q(:,:)
  integer,allocatable,save::old_classes(:,:),old_ids(:,:),old_sol(:,:),old_noh(:)
  integer(int2),allocatable,save::old_reverse(:,:)
  logical,save::boundary_suspended=.false.,boundary_pending=.false.
  real(wip),allocatable,target,save::boundary_translation(:,:,:),boundary_inverse(:,:)
  integer,allocatable,save::boundary_solute(:)
  integer(int2),allocatable,target,save::boundary_reverse(:,:)
  interface
    integer(c_int) function outgoing_mode() bind(C,name='genesis_atom_outgoing_mode')
      import c_int
    end function
    subroutine outgoing_c(view,data,box,origin,cell,grid,map,verify) bind(C,name='genesis_atom_outgoing')
      import particle_view,atom_records,c_double,c_int,c_ptr
      type(particle_view),intent(in)::view
      type(atom_records),intent(in)::data
      real(c_double),intent(in)::box(3),origin(3),cell(3)
      integer(c_int),intent(in)::grid(3)
      type(c_ptr),value::map
      integer(c_int),value::verify
    end subroutine
    subroutine outgoing_check_c() bind(C,name='genesis_atom_outgoing_check')
    end subroutine
    integer(c_int) function compare_c(a,b,n) bind(C,name="memcmp")
      import c_ptr,c_size_t,c_int
      type(c_ptr),value::a,b
      integer(c_size_t),value::n
    end function
    subroutine boundary_c(view,inverse,translation,box,origin,sites,epsilon) bind(C,name="genesis_incoming_boundary")
      import particle_view,c_ptr,c_double,c_int
      type(particle_view),intent(in)::view
      type(c_ptr),value::inverse,translation
      real(c_double),intent(in)::box(3),origin(3)
      integer(c_int),value::sites
      real(c_double),value::epsilon
    end subroutine
    integer(c_int) function mode_c() bind(C,name="genesis_incoming_mode")
      import c_int
    end function
    subroutine atom_c(view,data) bind(C,name="genesis_incoming_atoms")
      import particle_view,atom_records
      type(particle_view),intent(in)::view
      type(atom_records),intent(in)::data
    end subroutine
    subroutine water_c(view,data,sites,q,m,classes) bind(C,name="genesis_incoming_waters")
      import particle_view,records,c_int,c_ptr
      type(particle_view),intent(in)::view
      type(records),intent(in)::data
      integer(c_int),value::sites
      type(c_ptr),value::q,m,classes
    end subroutine
    subroutine group_c(view,data) bind(C,name="genesis_incoming_groups")
      import particle_view,records
      type(particle_view),intent(in)::view
      type(records),intent(in)::data
    end subroutine
  end interface
contains
  logical function incoming_gpu_enabled()
    incoming_gpu_enabled=mode_c()/=0
  end function
  logical function incoming_gpu_verify()
    incoming_gpu_verify=mode_c()==2.or.incoming_verify()
  end function
  function view_of(d,c) result(view)
    type(s_domain),target,intent(inout)::d
    type(s_constraints),target,optional,intent(in)::c
    type(particle_view)::view
    integer::width
    if(nproc_country/=1.or.d%num_cell_boundary/=0)error stop 'GPU incoming requires one rank without halos'
    if(wp/=c_double.or.wip/=c_double.or.int2/=c_int16_t)error stop 'GPU incoming requires double and 16-bit reverse map'
    width=size(d%coord,2)
    if(size(d%coord,1)/=3.or.size(d%velocity,2)/=width.or.size(d%charge,1)/=width.or. &
       size(d%mass,1)/=width.or.size(d%atom_cls_no,1)/=width.or.size(d%id_l2g,1)/=width.or. &
       size(d%id_l2g_solute,1)/=width.or.size(d%id_g2l,1)/=2)error stop 'GPU incoming pitch mismatch'
    if(width>32767.or.d%num_cell_local>32767.or.d%num_cell_local<1)error stop 'GPU incoming reverse-map width overflow'
    view%x=c_loc(d%coord);view%v=c_loc(d%velocity);view%q=c_loc(d%charge);view%m=c_loc(d%mass)
    view%classes=c_loc(d%atom_cls_no);view%ids=c_loc(d%id_l2g);view%solute=c_loc(d%id_l2g_solute)
    view%reverse=c_loc(d%id_g2l);view%ns=c_loc(d%num_solute);view%nw=c_loc(d%num_water);view%na=c_loc(d%num_atom)
    view%water=c_loc(d%water_list);view%width=width;view%cells=d%num_cell_local
    view%water_height=size(d%water_list,1);view%water_capacity=size(d%water_list,2)
    view%reverse_capacity=size(d%id_g2l,2)
    view%noh=c_null_ptr;view%groups=c_null_ptr;view%distance=c_null_ptr;view%hlocal=c_null_ptr
    view%connect=0;view%group_height=0;view%group_capacity=0;view%dist_height=0;view%dist_capacity=0
    if(present(c))then
      view%noh=c_loc(c%No_HGr);view%groups=c_loc(c%HGr_bond_list);view%distance=c_loc(c%HGr_bond_dist)
      view%hlocal=c_loc(c%HGr_local);view%connect=c%connect
      view%group_height=size(c%HGr_bond_list,1);view%group_capacity=size(c%HGr_bond_list,2)
      view%dist_height=size(c%HGr_bond_dist,1);view%dist_capacity=size(c%HGr_bond_dist,2)
      if(view%group_height<c%connect+1.or.view%dist_height<c%connect+1)error stop 'GPU incoming group height mismatch'
    end if
  end function
  subroutine incoming_water_gpu(sites,d)
    integer,intent(in)::sites
    type(s_domain),target,intent(inout)::d
    type(particle_view)::view
    type(records)::data
    view=view_of(d)
    if(sites<1.or.sites>view%water_height)error stop 'GPU incoming water height mismatch'
    data%move_real=c_loc(d%water%move_real);data%stay_real=c_loc(d%water%stay_real)
    data%move_int=c_loc(d%water%move_integer);data%stay_int=c_loc(d%water%stay_integer)
    data%move_count=c_loc(d%water%move);data%stay_count=c_loc(d%water%stay)
    data%move_real_rows=size(d%water%move_real,1);data%move_real_capacity=size(d%water%move_real,2)
    data%stay_real_rows=size(d%water%stay_real,1);data%stay_real_capacity=size(d%water%stay_real,2)
    data%move_int_rows=size(d%water%move_integer,1);data%move_int_capacity=size(d%water%move_integer,2)
    data%stay_int_rows=size(d%water%stay_integer,1);data%stay_int_capacity=size(d%water%stay_integer,2)
    if(min(data%move_real_rows,data%stay_real_rows)<6*sites.or. &
       min(data%move_int_rows,data%stay_int_rows)<sites)error stop 'GPU incoming water payload shape mismatch'
    call water_c(view,data,sites,c_loc(d%water%charge),c_loc(d%water%mass),c_loc(d%water%atom_cls_no))
  end subroutine
  subroutine incoming_group_gpu(c,d)
    type(s_constraints),target,intent(inout)::c
    type(s_domain),target,intent(inout)::d
    type(particle_view)::view
    type(records)::data
    view=view_of(d,c)
    data%move_real=c_loc(c%HGr_move_real);data%stay_real=c_loc(c%HGr_stay_real)
    data%move_int=c_loc(c%HGr_move_int);data%stay_int=c_loc(c%HGr_stay_int)
    data%move_count=c_loc(c%HGr_move);data%stay_count=c_loc(c%HGr_stay)
    data%move_real_rows=size(c%HGr_move_real,1);data%move_real_capacity=size(c%HGr_move_real,2)
    data%stay_real_rows=size(c%HGr_stay_real,1);data%stay_real_capacity=size(c%HGr_stay_real,2)
    data%move_int_rows=size(c%HGr_move_int,1);data%move_int_capacity=size(c%HGr_move_int,2)
    data%stay_int_rows=size(c%HGr_stay_int,1);data%stay_int_capacity=size(c%HGr_stay_int,2)
    if(min(data%move_real_rows,data%stay_real_rows)<9*(c%connect+1).or. &
       min(data%move_int_rows,data%stay_int_rows)<3*(c%connect+1))error stop 'GPU incoming group payload shape mismatch'
    call group_c(view,data)
  end subroutine
  subroutine incoming_atom_gpu(c,d)
    type(s_constraints),target,intent(inout)::c
    type(s_domain),target,intent(inout)::d
    type(particle_view)::view
    type(atom_records)::data
    view=view_of(d,c)
    data%reals=c_loc(d%buf_real);data%integers=c_loc(d%buf_integer)
    data%adds=c_loc(d%ptl_add);data%exits=c_loc(d%ptl_exit);data%indices=c_loc(d%ptl_exit_index)
    data%real_rows=size(d%buf_real,1);data%real_capacity=size(d%buf_real,2)
    data%int_rows=size(d%buf_integer,1);data%int_capacity=size(d%buf_integer,2)
    data%index_capacity=size(d%ptl_exit_index,1)
    if(data%real_rows<8.or.data%int_rows<3)error stop 'GPU incoming atom payload shape mismatch'
    call atom_c(view,data)
  end subroutine
  logical function outgoing_atom_verify()
    outgoing_atom_verify=outgoing_mode()==2.or.incoming_verify()
  end function
  subroutine outgoing_atom_check()
    call outgoing_check_c()
  end subroutine
  logical function outgoing_atom_gpu(b,c,d) result(handled)
    type(s_boundary),target,intent(in)::b
    type(s_constraints),target,intent(in)::c
    type(s_domain),target,intent(inout)::d
    type(particle_view)::view
    type(atom_records)::data
    handled=.false.
#ifdef USE_GPU
    if(outgoing_mode()==0)return
    if(.not.incoming_gpu_enabled())error stop 'GPU outgoing atoms require GPU incoming'
    view=view_of(d,c)
    data%reals=c_loc(d%buf_real);data%integers=c_loc(d%buf_integer)
    data%adds=c_loc(d%ptl_add);data%exits=c_loc(d%ptl_exit);data%indices=c_loc(d%ptl_exit_index)
    data%real_rows=size(d%buf_real,1);data%real_capacity=size(d%buf_real,2)
    data%int_rows=size(d%buf_integer,1);data%int_capacity=size(d%buf_integer,2)
    data%index_capacity=size(d%ptl_exit_index,1)
    call outgoing_c(view,data,[b%box_size_x,b%box_size_y,b%box_size_z], &
      [b%origin_x,b%origin_y,b%origin_z],[b%cell_size_x,b%cell_size_y,b%cell_size_z], &
      [b%num_cells_x,b%num_cells_y,b%num_cells_z],c_loc(d%cell_g2l),merge(1,0,outgoing_atom_verify()))
    handled=.true.
#endif
  end function
  subroutine incoming_atom_backup(d,c)
    type(s_domain),intent(in)::d
    type(s_constraints),intent(in)::c
    old_x=d%coord;old_v=d%velocity;old_m=d%mass;old_q=d%charge
    old_classes=d%atom_cls_no;old_ids=d%id_l2g;old_sol=d%id_l2g_solute;old_reverse=d%id_g2l;old_noh=c%No_HGr
  end subroutine
  subroutine incoming_atom_restore(d,c)
    type(s_domain),intent(inout)::d
    type(s_constraints),intent(inout)::c
    d%coord=old_x;d%velocity=old_v;d%mass=old_m;d%charge=old_q
    d%atom_cls_no=old_classes;d%id_l2g=old_ids;d%id_l2g_solute=old_sol;d%id_g2l=old_reverse;c%No_HGr=old_noh
  end subroutine
  subroutine incoming_boundary_suspend(value)
    logical,intent(in)::value
    boundary_suspended=value
  end subroutine
  logical function incoming_boundary_gpu(d,c,b) result(handled)
    type(s_domain),target,intent(inout)::d
    type(s_constraints),target,intent(inout)::c
    type(s_boundary),intent(in)::b
    type(particle_view)::view
    integer::sites
    real(c_double)::box(3),origin(3)
    handled=.false.;boundary_pending=.false.
    if(boundary_suspended.or..not.incoming_gpu_enabled())return
    view=view_of(d,c)
    select case(c%water_type)
    case(TIP1);sites=1
    case(TIP2);sites=2
    case(TIP3);sites=3
    case(TIP4);sites=4
    case default;error stop 'GPU bookkeeping unsupported water type'
    end select
    if(size(d%inv_mass,1)/=view%width.or.size(d%trans_vec,2)/=view%width)error stop 'GPU bookkeeping pitch mismatch'
    box=[b%box_size_x,b%box_size_y,b%box_size_z];origin=[b%origin_x,b%origin_y,b%origin_z]
    call boundary_c(view,c_loc(d%inv_mass),c_loc(d%trans_vec),box,origin,sites,real(EPS,c_double))
    handled=.not.incoming_gpu_verify()
    if(.not.handled)then
      boundary_translation=d%trans_vec;boundary_inverse=d%inv_mass
      boundary_solute=d%num_solute;boundary_reverse=d%id_g2l;boundary_pending=.true.
    end if
  end function
  subroutine incoming_boundary_check(d)
    type(s_domain),target,intent(in)::d
    integer::cell,n,atom,id
    if(.not.boundary_pending)return
    do cell=1,d%num_cell_local
      n=d%num_atom(cell)
      if(boundary_solute(cell)/=d%num_solute(cell))error stop 'GPU bookkeeping solute count mismatch'
      if(n>0)then
        if(compare_c(c_loc(boundary_translation(1,1,cell)),c_loc(d%trans_vec(1,1,cell)),int(24*n,c_size_t))/=0) &
          error stop 'GPU bookkeeping translations differ from CPU'
        if(compare_c(c_loc(boundary_inverse(1,cell)),c_loc(d%inv_mass(1,cell)),int(8*n,c_size_t))/=0) &
          error stop 'GPU bookkeeping inverse mass differs from CPU'
      end if
      do atom=1,d%num_solute(cell)
        id=d%id_l2g_solute(atom,cell)
        if(compare_c(c_loc(boundary_reverse(1,id)),c_loc(d%id_g2l(1,id)),4_c_size_t)/=0) &
          error stop 'GPU bookkeeping reverse map differs from CPU'
      end do
    end do
    boundary_pending=.false.
    write(*,'(a)') 'GENESIS_INCOMING_GPU: exact cell bookkeeping pass'
  end subroutine
end module
