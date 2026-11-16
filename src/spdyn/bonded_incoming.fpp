module bonded_incoming_mod
  use iso_c_binding
  use constants_mod
  use mpi_parallel_mod
  use sp_domain_str_mod
  use sp_enefunc_str_mod
  use bonded_migration_mod,only:bonded_verify,bonded_identity_fields
  implicit none
  private
  public::bonded_improper_bootstrap
  public::bonded_incoming_begin,bonded_incoming_finish,bonded_incoming_reference
  type,bind(C)::move_field
    type(c_ptr)::target,source
    integer(c_int)::rows,capacity,source_rows,source_capacity,source_offset,bytes
  end type
  type,bind(C)::move_view
    type(c_ptr)::counts,adds,exits,indices
    integer(c_int)::cells,index_capacity,fields
    type(move_field)::field(12)
  end type
  logical,save::reference_needed=.false.
  interface
    integer(c_int) function hbm_ready_c(kind) bind(C,name='genesis_bonded_hbm_ready')
      import c_int
      integer(c_int),value::kind
    end function
    subroutine initialize_c(view,kind) bind(C,name='genesis_bonded_hbm_initialize')
      import move_view,c_int
      type(move_view),intent(in)::view
      integer(c_int),value::kind
    end subroutine
    integer(c_int) function mode_c() bind(C,name="genesis_bonded_incoming_mode")
      import c_int
    end function
    subroutine start_c(view,kind,verify) bind(C,name="genesis_bonded_incoming_start")
      import move_view,c_int
      type(move_view),intent(in)::view
      integer(c_int),value::kind,verify
    end subroutine
    subroutine finish_c(kind) bind(C,name="genesis_bonded_incoming_finish")
      import c_int
      integer(c_int),value::kind
    end subroutine
  end interface
contains
  logical function bonded_incoming_reference()
    bonded_incoming_reference=reference_needed
  end function
  subroutine bonded_incoming_finish(kind)
    integer,intent(in)::kind
    call finish_c(kind)
    reference_needed=.false.
  end subroutine
  subroutine add(view,target,rows,capacity,source,source_rows,source_capacity,offset,bytes)
    type(move_view),intent(inout)::view
    type(c_ptr),intent(in)::target,source
    integer,intent(in)::rows,capacity,source_rows,source_capacity,offset,bytes
    if(rows<=0.or.capacity<=0.or.source_rows<offset+rows.or.source_capacity<=0)error stop 'Invalid bonded field shape'
    view%fields=view%fields+1
    if(view%fields>12)error stop 'Too many bonded record fields'
    view%field(view%fields)=move_field(target,source,rows,capacity,source_rows,source_capacity,offset,bytes)
  end subroutine
  subroutine bonded_improper_bootstrap(d,e)
    type(s_domain),intent(in)::d
    type(s_enefunc),target,intent(inout)::e
    logical::prepared
    if(hbm_ready_c(4)==0)prepared=bonded_incoming_begin(4,d,e,.true.)
  end subroutine
  logical function bonded_incoming_begin(kind,d,e,initialize_only) result(active)
    integer,intent(in)::kind
    type(s_domain),intent(in)::d
    type(s_enefunc),target,intent(inout)::e
    logical,optional,intent(in)::initialize_only
    type(move_view)::view
    type(c_ptr)::identity,incoming
    integer::capacity,incoming_capacity
    active=.false.;reference_needed=.false.
    if(mode_c()==0)return
    if(nproc_country/=1.or.d%num_cell_boundary/=0.or.e%forcefield/=ForcefieldCHARMM) &
      error stop 'GPU bonded incoming requires one-rank CHARMM without halos'
    active=.true.
    if(kind==5.and.e%num_rb_dihe_all==0)return
    if(kind==4.and.e%num_impr_all==0)return
    view%cells=d%num_cell_local;view%fields=0
    view%field=move_field(c_null_ptr,c_null_ptr,0,0,0,0,0,0)
    select case(kind)
    case(1)
      view%counts=c_loc(e%num_bond);view%adds=c_loc(e%bond_add);view%exits=c_loc(e%bond_exit)
      view%indices=c_loc(e%bond_exit_index);view%index_capacity=size(e%bond_exit_index,1)
      call add(view,c_loc(e%bond_list),size(e%bond_list,1),size(e%bond_list,2), &
        c_loc(e%buf_bond_integer),size(e%buf_bond_integer,1),size(e%buf_bond_integer,2),0,4)
      call add(view,c_loc(e%bond_kind),1,size(e%bond_kind,1), &
        c_loc(e%buf_bond_integer),size(e%buf_bond_integer,1),size(e%buf_bond_integer,2),2,1)
      call add(view,c_loc(e%bond_pbc),1,size(e%bond_pbc,1), &
        c_loc(e%buf_bond_integer),size(e%buf_bond_integer,1),size(e%buf_bond_integer,2),3,4)
      call add(view,c_loc(e%bond_dist_min),1,size(e%bond_dist_min,1), &
        c_loc(e%buf_bond_real),size(e%buf_bond_real,1),size(e%buf_bond_real,2),0,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%bond_force_const),1,size(e%bond_force_const,1), &
        c_loc(e%buf_bond_real),size(e%buf_bond_real,1),size(e%buf_bond_real,2),1,storage_size(1.0_wp)/8)
    case(2)
      view%counts=c_loc(e%num_angle);view%adds=c_loc(e%angle_add);view%exits=c_loc(e%angle_exit)
      view%indices=c_loc(e%angle_exit_index);view%index_capacity=size(e%angle_exit_index,1)
      call add(view,c_loc(e%angle_list),size(e%angle_list,1),size(e%angle_list,2), &
        c_loc(e%buf_angle_integer),size(e%buf_angle_integer,1),size(e%buf_angle_integer,2),0,4)
      call add(view,c_loc(e%angle_kind),1,size(e%angle_kind,1), &
        c_loc(e%buf_angle_integer),size(e%buf_angle_integer,1),size(e%buf_angle_integer,2),3,1)
      call add(view,c_loc(e%angle_pbc),size(e%angle_pbc,1),size(e%angle_pbc,2), &
        c_loc(e%buf_angle_integer),size(e%buf_angle_integer,1),size(e%buf_angle_integer,2),4,4)
      call add(view,c_loc(e%angle_theta_min),1,size(e%angle_theta_min,1), &
        c_loc(e%buf_angle_real),size(e%buf_angle_real,1),size(e%buf_angle_real,2),0,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%angle_force_const),1,size(e%angle_force_const,1), &
        c_loc(e%buf_angle_real),size(e%buf_angle_real,1),size(e%buf_angle_real,2),1,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%urey_rmin),1,size(e%urey_rmin,1), &
        c_loc(e%buf_angle_real),size(e%buf_angle_real,1),size(e%buf_angle_real,2),2,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%urey_force_const),1,size(e%urey_force_const,1), &
        c_loc(e%buf_angle_real),size(e%buf_angle_real,1),size(e%buf_angle_real,2),3,storage_size(1.0_wp)/8)
    case(3)
      view%counts=c_loc(e%num_dihedral);view%adds=c_loc(e%dihed_add);view%exits=c_loc(e%dihed_exit)
      view%indices=c_loc(e%dihed_exit_index);view%index_capacity=size(e%dihed_exit_index,1)
      call add(view,c_loc(e%dihe_list),size(e%dihe_list,1),size(e%dihe_list,2), &
        c_loc(e%buf_dihed_integer),size(e%buf_dihed_integer,1),size(e%buf_dihed_integer,2),0,4)
      call add(view,c_loc(e%dihe_periodicity),1,size(e%dihe_periodicity,1), &
        c_loc(e%buf_dihed_integer),size(e%buf_dihed_integer,1),size(e%buf_dihed_integer,2),4,4)
      call add(view,c_loc(e%dihe_kind),1,size(e%dihe_kind,1), &
        c_loc(e%buf_dihed_integer),size(e%buf_dihed_integer,1),size(e%buf_dihed_integer,2),5,1)
      call add(view,c_loc(e%dihe_pbc),size(e%dihe_pbc,1),size(e%dihe_pbc,2), &
        c_loc(e%buf_dihed_integer),size(e%buf_dihed_integer,1),size(e%buf_dihed_integer,2),6,4)
      call add(view,c_loc(e%dihe_phase),1,size(e%dihe_phase,1), &
        c_loc(e%buf_dihed_real),size(e%buf_dihed_real,1),size(e%buf_dihed_real,2),0,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%dihe_force_const),1,size(e%dihe_force_const,1), &
        c_loc(e%buf_dihed_real),size(e%buf_dihed_real,1),size(e%buf_dihed_real,2),1,storage_size(1.0_wp)/8)
    case(4)
      view%counts=c_loc(e%num_improper);view%adds=c_loc(e%impr_add);view%exits=c_loc(e%impr_exit)
      view%indices=c_loc(e%impr_exit_index);view%index_capacity=size(e%impr_exit_index,1)
      call add(view,c_loc(e%impr_list),size(e%impr_list,1),size(e%impr_list,2), &
        c_loc(e%buf_impr_integer),size(e%buf_impr_integer,1),size(e%buf_impr_integer,2),0,4)
      call add(view,c_loc(e%impr_periodicity),1,size(e%impr_periodicity,1), &
        c_loc(e%buf_impr_integer),size(e%buf_impr_integer,1),size(e%buf_impr_integer,2),4,4)
      call add(view,c_loc(e%impr_pbc),size(e%impr_pbc,1),size(e%impr_pbc,2), &
        c_loc(e%buf_impr_integer),size(e%buf_impr_integer,1),size(e%buf_impr_integer,2),5,4)
      call add(view,c_loc(e%impr_phase),1,size(e%impr_phase,1), &
        c_loc(e%buf_impr_real),size(e%buf_impr_real,1),size(e%buf_impr_real,2),0,storage_size(1.0_wp)/8)
      call add(view,c_loc(e%impr_force_const),1,size(e%impr_force_const,1), &
        c_loc(e%buf_impr_real),size(e%buf_impr_real,1),size(e%buf_impr_real,2),1,storage_size(1.0_wp)/8)
    case(5)
      view%counts=c_loc(e%num_rb_dihedral);view%adds=c_loc(e%rb_dihed_add);view%exits=c_loc(e%rb_dihed_exit)
      view%indices=c_loc(e%rb_dihed_exit_index);view%index_capacity=size(e%rb_dihed_exit_index,1)
      call add(view,c_loc(e%rb_dihe_list),size(e%rb_dihe_list,1),size(e%rb_dihe_list,2), &
        c_loc(e%buf_rb_dihed_integer),size(e%buf_rb_dihed_integer,1),size(e%buf_rb_dihed_integer,2),0,4)
      call add(view,c_loc(e%rb_dihe_pbc),size(e%rb_dihe_pbc,1),size(e%rb_dihe_pbc,2), &
        c_loc(e%buf_rb_dihed_integer),size(e%buf_rb_dihed_integer,1),size(e%buf_rb_dihed_integer,2),4,4)
      call add(view,c_loc(e%rb_dihe_c),size(e%rb_dihe_c,1),size(e%rb_dihe_c,2), &
        c_loc(e%buf_rb_dihed_real),size(e%buf_rb_dihed_real,1),size(e%buf_rb_dihed_real,2),0,storage_size(1.0_wp)/8)
    case default
      error stop 'Unknown bonded incoming kind'
    end select
    if(kind<=3)then
      call bonded_identity_fields(kind,identity,incoming,capacity,incoming_capacity)
      call add(view,identity,1,capacity,incoming,1,incoming_capacity,0,4)
    end if
    if(c_sizeof(view)/=528_c_size_t)error stop 'Bonded incoming C layout mismatch'
    if(present(initialize_only))then
      if(initialize_only)then
        call initialize_c(view,kind)
        return
      end if
    end if
    reference_needed=mode_c()==2.or.bonded_verify()
    call start_c(view,kind,merge(1,0,reference_needed))
  end function
end module
