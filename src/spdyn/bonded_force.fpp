#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module bonded_force_mod
  use cudafor
  use iso_c_binding
  use constants_mod
  use mpi_parallel_mod
  use sp_domain_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod,only:VDWPME
  use sp_pairlist_str_mod
  use coordinate_pack_mod,only:coordinate_copy_c
  use genesis_step_state_mod,only:step_export_force_coord
  use force_assembly_mod,only:force_assembly_enabled
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public :: bonded_enabled,bonded_cpu_reference,bonded_prepare,bonded_finish
  integer::cells=0,width=0,atoms=0
  integer,device::g_atoms,g_notation
  real(8),device,allocatable::g_box(:)
  real(8),device,allocatable::g_coord(:,:,:),g_force(:,:,:),g_pbc(:,:,:)
  integer(int2),device,allocatable::g_id(:,:)
  integer,device,allocatable::g_class(:)
  real(8),allocatable,pinned::host_parts(:,:,:)
  real(8),allocatable::before_values(:,:)
  integer(kind=cuda_stream_kind)::producer_stream
  real(8),allocatable,pinned::snapshot_coord(:,:,:),snapshot_pbc(:,:,:),snapshot_box(:)
  logical::first=.true.
  interface
    subroutine register_force(values) bind(C,name='genesis_bonded_register')
      import c_double
      real(c_double),device,intent(in)::values(*)
    end subroutine
  end interface
  integer,device,allocatable::d_bond_bond_pbc(:)
  integer,device,allocatable::d_bond_bondlist(:,:)
  real(wp),device,allocatable::d_bond_fc(:)
  integer,device,allocatable::d_bond_nbond(:)
  real(wp),device,allocatable::d_bond_r0(:)
  real(8),device,allocatable::out_bond(:,:),sum_bond(:,:)
  integer,device,allocatable::d_angle_angl_pbc(:,:)
  integer,device,allocatable::d_angle_anglelist(:,:)
  real(wp),device,allocatable::d_angle_fc(:)
  real(wp),device,allocatable::d_angle_fc_ub(:)
  integer,device,allocatable::d_angle_nangle(:)
  real(wp),device,allocatable::d_angle_r0_ub(:)
  real(wp),device,allocatable::d_angle_theta0(:)
  real(8),device,allocatable::out_angle(:,:),sum_angle(:,:)
  integer,device,allocatable::d_dihed_dihe_pbc(:,:)
  integer,device,allocatable::d_dihed_dihelist(:,:)
  real(wp),device,allocatable::d_dihed_fc(:)
  integer,device,allocatable::d_dihed_ndihe(:)
  integer,device,allocatable::d_dihed_nperiod(:)
  real(wp),device,allocatable::d_dihed_phase(:)
  real(8),device,allocatable::out_dihed(:,:),sum_dihed(:,:)
  real(wp),device,allocatable::d_improper_fc(:)
  integer,device,allocatable::d_improper_impr_pbc(:,:)
  integer,device,allocatable::d_improper_imprlist(:,:)
  integer,device,allocatable::d_improper_nimp(:)
  real(wp),device,allocatable::d_improper_phase(:)
  real(8),device,allocatable::out_improper(:,:),sum_improper(:,:)
  integer(1),device,allocatable::d_nb14_cell_move(:,:)
  real(wp),device::d_nb14_cutoff
  real(wp),device::d_nb14_density
  integer,device,allocatable::d_nb14_nb14_calc_list(:,:)
  real(wp),device,allocatable::d_nb14_nb14_lj12(:,:)
  real(wp),device,allocatable::d_nb14_nb14_lj6(:,:)
  integer,device,allocatable::d_nb14_num_nb14_calc(:)
  integer,device,allocatable::d_nb14_start_atom(:)
  real(wp),device,allocatable::d_nb14_table_ene(:)
  real(wp),device,allocatable::d_nb14_table_grad(:)
  real(8),device,allocatable::out_nb14(:,:),sum_nb14(:,:)
  integer(1),device,allocatable::d_corr_cell_move(:,:)
  real(wp),device::d_corr_cutoff
  real(wp),device::d_corr_density
  integer,device,allocatable::d_corr_nonb_excl_list(:,:)
  integer,device,allocatable::d_corr_num_nonb_excl(:)
  integer,device,allocatable::d_corr_start_atom(:)
  real(wp),device,allocatable::d_corr_table_decor(:)
  real(wp),device,allocatable::d_corr_table_ecor(:)
  real(8),device,allocatable::out_corr(:,:),sum_corr(:,:)

  integer::capacity_bond=0
  integer,device::live_bond
  integer,device,allocatable::offset_bond(:),owner_bond(:),local_bond(:)
  integer,allocatable,pinned::h_offset_bond(:),h_owner_bond(:),h_local_bond(:)
  integer,allocatable,pinned::h_bond_bondlist(:,:)
  real(wp),allocatable,pinned::h_bond_fc(:)
  real(wp),allocatable,pinned::h_bond_r0(:)
  integer,allocatable,pinned::h_bond_bond_pbc(:)
  integer::capacity_angle=0
  integer,device::live_angle
  integer,device,allocatable::offset_angle(:),owner_angle(:),local_angle(:)
  integer,allocatable,pinned::h_offset_angle(:),h_owner_angle(:),h_local_angle(:)
  integer,allocatable,pinned::h_angle_anglelist(:,:)
  real(wp),allocatable,pinned::h_angle_fc(:)
  real(wp),allocatable,pinned::h_angle_theta0(:)
  real(wp),allocatable,pinned::h_angle_fc_ub(:)
  real(wp),allocatable,pinned::h_angle_r0_ub(:)
  integer,allocatable,pinned::h_angle_angl_pbc(:,:)
  integer::capacity_dihed=0
  integer,device::live_dihed
  integer,device,allocatable::offset_dihed(:),owner_dihed(:),local_dihed(:)
  integer,allocatable,pinned::h_offset_dihed(:),h_owner_dihed(:),h_local_dihed(:)
  integer,allocatable,pinned::h_dihed_dihelist(:,:)
  real(wp),allocatable,pinned::h_dihed_fc(:)
  integer,allocatable,pinned::h_dihed_nperiod(:)
  real(wp),allocatable,pinned::h_dihed_phase(:)
  integer,allocatable,pinned::h_dihed_dihe_pbc(:,:)
  integer::capacity_improper=0
  integer,device::live_improper
  integer,device,allocatable::offset_improper(:),owner_improper(:),local_improper(:)
  integer,allocatable,pinned::h_offset_improper(:),h_owner_improper(:),h_local_improper(:)
  integer,allocatable,pinned::h_improper_imprlist(:,:)
  real(wp),allocatable,pinned::h_improper_fc(:)
  real(wp),allocatable,pinned::h_improper_phase(:)
  integer,allocatable,pinned::h_improper_impr_pbc(:,:)
  integer::capacity_nb14=0
  integer,device::live_nb14
  integer,device,allocatable::offset_nb14(:),owner_nb14(:),local_nb14(:)
  integer,allocatable,pinned::h_offset_nb14(:),h_owner_nb14(:),h_local_nb14(:)
  integer(1),allocatable,pinned::h_nb14_cell_move(:,:)
  integer,allocatable,pinned::h_nb14_nb14_calc_list(:,:)
  integer::capacity_corr=0
  integer,device::live_corr
  integer,device,allocatable::offset_corr(:),owner_corr(:),local_corr(:)
  integer,allocatable,pinned::h_offset_corr(:),h_owner_corr(:),h_local_corr(:)
  integer(1),allocatable,pinned::h_corr_cell_move(:,:)
  integer,allocatable,pinned::h_corr_nonb_excl_list(:,:)
contains
  logical function bonded_enabled()
    integer,save::selected=-1
    integer::status,code
    character(16)::value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_FORCE_BONDED',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0.or.selected<0.or.selected>1) error stop 'GENESIS_FORCE_BONDED must be 0 or 1'
      end if
    end if
    bonded_enabled=selected==1
    if(bonded_enabled) then
      if(.not.force_assembly_enabled()) error stop 'GPU bonded terms require GPU force assembly'
    end if
  end function

  logical function bonded_verify()
    integer::status
    character(16)::value
    call get_environment_variable('GENESIS_FORCE_VERIFY',value,status=status)
    bonded_verify=status==0.and.trim(value)=='1'
  end function

  logical function bonded_cpu_reference()
    bonded_cpu_reference=.not.bonded_enabled().or.bonded_verify()
  end function

  subroutine checked(status)
    integer,intent(in)::status
    if(status/=cudaSuccess) then
      print *, 'GENESIS bonded CUDA error ',status
      error stop 'GPU bonded runtime error'
    end if
  end subroutine

  subroutine values_pack(values,vir,b,a,u,d,i,e,v)
    real(8),intent(out)::values(:,:)
    real(dp),intent(in)::vir(:,:,:),b(:),a(:),u(:),d(:),i(:),e(:),v(:)
    values(1,:)=b;values(2,:)=a;values(3,:)=u;values(4,:)=d;values(5,:)=i
    values(6,:)=e;values(7,:)=v
    values(8,:)=vir(1,1,:);values(9,:)=vir(2,2,:);values(10,:)=vir(3,3,:)
  end subroutine

  subroutine bonded_prepare(domain,enefunc,pairlist,coord,pbc,classes,vir,b,a,u,d,i,e,v)
    type(s_domain),intent(in)::domain
    type(s_enefunc),intent(in)::enefunc
    type(s_pairlist),intent(in)::pairlist
    real(wip),intent(in)::coord(:,:,:)
    real(wp),intent(in)::pbc(:,:,:)
    integer,intent(in)::classes(:)
    real(dp),intent(in)::vir(:,:,:),b(:),a(:),u(:),d(:),i(:),e(:),v(:)
    integer::icell,entry,slot,total
    integer(kind=cuda_stream_kind)::previous_stream
    logical::copied
    if(.not.real_calc.or.enefunc%vdw==VDWPME.or.enefunc%nonb_limiter.or.enefunc%local_restraint.or.enefunc%restraint.or.enefunc%use_efield) &
      error stop 'GPU bonded path requires CHARMM without LJ-PME, limiters, restraints or fields'
    if(domain%num_cell_boundary/=0) error stop 'GPU bonded halo unsupported'
    if(first) then
      call checked(cudaStreamCreateWithFlags(producer_stream,cudaStreamNonBlocking))
      allocate(g_box(3))
      cells=domain%num_cell_local;width=size(coord,2);atoms=domain%num_atom_domain
      capacity_bond=max(1,sum(enefunc%num_bond(1:cells)))
      allocate(offset_bond(cells),owner_bond(capacity_bond),local_bond(capacity_bond))
      allocate(h_offset_bond(cells),h_owner_bond(capacity_bond),h_local_bond(capacity_bond))
      allocate(h_bond_bondlist(size(enefunc%bond_list,1),capacity_bond))
      allocate(h_bond_fc(capacity_bond))
      allocate(h_bond_r0(capacity_bond))
      allocate(h_bond_bond_pbc(capacity_bond))
      capacity_angle=max(1,sum(enefunc%num_angle(1:cells)))
      allocate(offset_angle(cells),owner_angle(capacity_angle),local_angle(capacity_angle))
      allocate(h_offset_angle(cells),h_owner_angle(capacity_angle),h_local_angle(capacity_angle))
      allocate(h_angle_anglelist(size(enefunc%angle_list,1),capacity_angle))
      allocate(h_angle_fc(capacity_angle))
      allocate(h_angle_theta0(capacity_angle))
      allocate(h_angle_fc_ub(capacity_angle))
      allocate(h_angle_r0_ub(capacity_angle))
      allocate(h_angle_angl_pbc(size(enefunc%angle_pbc,1),capacity_angle))
      capacity_dihed=max(1,sum(enefunc%num_dihedral(1:cells)))
      allocate(offset_dihed(cells),owner_dihed(capacity_dihed),local_dihed(capacity_dihed))
      allocate(h_offset_dihed(cells),h_owner_dihed(capacity_dihed),h_local_dihed(capacity_dihed))
      allocate(h_dihed_dihelist(size(enefunc%dihe_list,1),capacity_dihed))
      allocate(h_dihed_fc(capacity_dihed))
      allocate(h_dihed_nperiod(capacity_dihed))
      allocate(h_dihed_phase(capacity_dihed))
      allocate(h_dihed_dihe_pbc(size(enefunc%dihe_pbc,1),capacity_dihed))
      capacity_improper=max(1,sum(enefunc%num_improper(1:cells)))
      allocate(offset_improper(cells),owner_improper(capacity_improper),local_improper(capacity_improper))
      allocate(h_offset_improper(cells),h_owner_improper(capacity_improper),h_local_improper(capacity_improper))
      allocate(h_improper_imprlist(size(enefunc%impr_list,1),capacity_improper))
      allocate(h_improper_fc(capacity_improper))
      allocate(h_improper_phase(capacity_improper))
      allocate(h_improper_impr_pbc(size(enefunc%impr_pbc,1),capacity_improper))
      capacity_nb14=max(1,sum(enefunc%num_nb14_calc(1:cells)))
      allocate(offset_nb14(cells),owner_nb14(capacity_nb14),local_nb14(capacity_nb14))
      allocate(h_offset_nb14(cells),h_owner_nb14(capacity_nb14),h_local_nb14(capacity_nb14))
      allocate(h_nb14_cell_move(size(domain%cell_move,1),capacity_nb14))
      allocate(h_nb14_nb14_calc_list(size(enefunc%nb14_calc_list,1),capacity_nb14))
      capacity_corr=max(1,sum(enefunc%num_nonb_excl(1:cells)))
      allocate(offset_corr(cells),owner_corr(capacity_corr),local_corr(capacity_corr))
      allocate(h_offset_corr(cells),h_owner_corr(capacity_corr),h_local_corr(capacity_corr))
      allocate(h_corr_cell_move(size(domain%cell_move,1),capacity_corr))
      allocate(h_corr_nonb_excl_list(size(enefunc%nonb_excl_list,1),capacity_corr))
      allocate(snapshot_coord(3,width,cells),snapshot_pbc(4*atoms,1,1),snapshot_box(3))
      allocate(g_coord(3,width,cells),g_force(3,width,cells),g_pbc(4*atoms,1,1))
      allocate(g_id(size(domain%id_g2l,1),size(domain%id_g2l,2)),g_class(atoms))
      allocate(host_parts(10,cells,6),before_values(10,nthread))
      allocate(d_bond_bond_pbc(capacity_bond))
      allocate(d_bond_bondlist(size(enefunc%bond_list,1),capacity_bond))
      allocate(d_bond_fc(capacity_bond))
      allocate(d_bond_nbond(size(enefunc%num_bond,1)))
      allocate(d_bond_r0(capacity_bond))
      allocate(out_bond(10,capacity_bond),sum_bond(10,cells))
      allocate(d_angle_angl_pbc(size(enefunc%angle_pbc,1),capacity_angle))
      allocate(d_angle_anglelist(size(enefunc%angle_list,1),capacity_angle))
      allocate(d_angle_fc(capacity_angle))
      allocate(d_angle_fc_ub(capacity_angle))
      allocate(d_angle_nangle(size(enefunc%num_angle,1)))
      allocate(d_angle_r0_ub(capacity_angle))
      allocate(d_angle_theta0(capacity_angle))
      allocate(out_angle(10,capacity_angle),sum_angle(10,cells))
      allocate(d_dihed_dihe_pbc(size(enefunc%dihe_pbc,1),capacity_dihed))
      allocate(d_dihed_dihelist(size(enefunc%dihe_list,1),capacity_dihed))
      allocate(d_dihed_fc(capacity_dihed))
      allocate(d_dihed_ndihe(size(enefunc%num_dihedral,1)))
      allocate(d_dihed_nperiod(capacity_dihed))
      allocate(d_dihed_phase(capacity_dihed))
      allocate(out_dihed(10,capacity_dihed),sum_dihed(10,cells))
      allocate(d_improper_fc(capacity_improper))
      allocate(d_improper_impr_pbc(size(enefunc%impr_pbc,1),capacity_improper))
      allocate(d_improper_imprlist(size(enefunc%impr_list,1),capacity_improper))
      allocate(d_improper_nimp(size(enefunc%num_improper,1)))
      allocate(d_improper_phase(capacity_improper))
      allocate(out_improper(10,capacity_improper),sum_improper(10,cells))
      allocate(d_nb14_cell_move(size(domain%cell_move,1),capacity_nb14))
      allocate(d_nb14_nb14_calc_list(size(enefunc%nb14_calc_list,1),capacity_nb14))
      allocate(d_nb14_nb14_lj12(size(enefunc%nb14_lj12,1),size(enefunc%nb14_lj12,2)))
      allocate(d_nb14_nb14_lj6(size(enefunc%nb14_lj6,1),size(enefunc%nb14_lj6,2)))
      allocate(d_nb14_num_nb14_calc(size(enefunc%num_nb14_calc,1)))
      allocate(d_nb14_start_atom(size(domain%start_atom,1)))
      allocate(d_nb14_table_ene(size(enefunc%table%table_ene,1)))
      allocate(d_nb14_table_grad(size(enefunc%table%table_grad,1)))
      allocate(out_nb14(10,capacity_nb14),sum_nb14(10,cells))
      allocate(d_corr_cell_move(size(domain%cell_move,1),capacity_corr))
      allocate(d_corr_nonb_excl_list(size(enefunc%nonb_excl_list,1),capacity_corr))
      allocate(d_corr_num_nonb_excl(size(enefunc%num_nonb_excl,1)))
      allocate(d_corr_start_atom(size(domain%start_atom,1)))
      allocate(d_corr_table_decor(size(enefunc%table%table_decor,1)))
      allocate(d_corr_table_ecor(size(enefunc%table%table_ecor,1)))
      allocate(out_corr(10,capacity_corr),sum_corr(10,cells))

    end if
    if(first) then
      call checked(cudaforSetDefaultStream(g_box,producer_stream))
      call checked(cudaforSetDefaultStream(g_coord,producer_stream))
      call checked(cudaforSetDefaultStream(g_force,producer_stream))
      call checked(cudaforSetDefaultStream(g_pbc,producer_stream))
      call checked(cudaforSetDefaultStream(g_id,producer_stream))
      call checked(cudaforSetDefaultStream(g_class,producer_stream))
      call checked(cudaforSetDefaultStream(d_bond_bond_pbc,producer_stream))
      call checked(cudaforSetDefaultStream(d_bond_bondlist,producer_stream))
      call checked(cudaforSetDefaultStream(d_bond_fc,producer_stream))
      call checked(cudaforSetDefaultStream(d_bond_nbond,producer_stream))
      call checked(cudaforSetDefaultStream(d_bond_r0,producer_stream))
      call checked(cudaforSetDefaultStream(out_bond,producer_stream))
      call checked(cudaforSetDefaultStream(sum_bond,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_angl_pbc,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_anglelist,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_fc,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_fc_ub,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_nangle,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_r0_ub,producer_stream))
      call checked(cudaforSetDefaultStream(d_angle_theta0,producer_stream))
      call checked(cudaforSetDefaultStream(out_angle,producer_stream))
      call checked(cudaforSetDefaultStream(sum_angle,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_dihe_pbc,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_dihelist,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_fc,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_ndihe,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_nperiod,producer_stream))
      call checked(cudaforSetDefaultStream(d_dihed_phase,producer_stream))
      call checked(cudaforSetDefaultStream(out_dihed,producer_stream))
      call checked(cudaforSetDefaultStream(sum_dihed,producer_stream))
      call checked(cudaforSetDefaultStream(d_improper_fc,producer_stream))
      call checked(cudaforSetDefaultStream(d_improper_impr_pbc,producer_stream))
      call checked(cudaforSetDefaultStream(d_improper_imprlist,producer_stream))
      call checked(cudaforSetDefaultStream(d_improper_nimp,producer_stream))
      call checked(cudaforSetDefaultStream(d_improper_phase,producer_stream))
      call checked(cudaforSetDefaultStream(out_improper,producer_stream))
      call checked(cudaforSetDefaultStream(sum_improper,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_cell_move,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_nb14_calc_list,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_nb14_lj12,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_nb14_lj6,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_num_nb14_calc,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_start_atom,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_table_ene,producer_stream))
      call checked(cudaforSetDefaultStream(d_nb14_table_grad,producer_stream))
      call checked(cudaforSetDefaultStream(out_nb14,producer_stream))
      call checked(cudaforSetDefaultStream(sum_nb14,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_cell_move,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_nonb_excl_list,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_num_nonb_excl,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_start_atom,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_table_decor,producer_stream))
      call checked(cudaforSetDefaultStream(d_corr_table_ecor,producer_stream))
      call checked(cudaforSetDefaultStream(out_corr,producer_stream))
      call checked(cudaforSetDefaultStream(sum_corr,producer_stream))
      call checked(cudaforSetDefaultStream(offset_bond,producer_stream))
      call checked(cudaforSetDefaultStream(owner_bond,producer_stream))
      call checked(cudaforSetDefaultStream(local_bond,producer_stream))
      call checked(cudaforSetDefaultStream(offset_angle,producer_stream))
      call checked(cudaforSetDefaultStream(owner_angle,producer_stream))
      call checked(cudaforSetDefaultStream(local_angle,producer_stream))
      call checked(cudaforSetDefaultStream(offset_dihed,producer_stream))
      call checked(cudaforSetDefaultStream(owner_dihed,producer_stream))
      call checked(cudaforSetDefaultStream(local_dihed,producer_stream))
      call checked(cudaforSetDefaultStream(offset_improper,producer_stream))
      call checked(cudaforSetDefaultStream(owner_improper,producer_stream))
      call checked(cudaforSetDefaultStream(local_improper,producer_stream))
      call checked(cudaforSetDefaultStream(offset_nb14,producer_stream))
      call checked(cudaforSetDefaultStream(owner_nb14,producer_stream))
      call checked(cudaforSetDefaultStream(local_nb14,producer_stream))
      call checked(cudaforSetDefaultStream(offset_corr,producer_stream))
      call checked(cudaforSetDefaultStream(owner_corr,producer_stream))
      call checked(cudaforSetDefaultStream(local_corr,producer_stream))
    end if
    if(cells/=domain%num_cell_local.or.width/=size(coord,2).or.atoms/=domain%num_atom_domain) &
      error stop 'GPU bonded shape changed'
    previous_stream=cudaforGetDefaultStream()
    call checked(cudaforSetDefaultStream(stream=producer_stream))
    if(first.or.pairlist%univ_update/=0) then
      g_id=domain%id_g2l;g_class=classes(1:atoms)
      g_atoms=atoms;g_notation=enefunc%notation_14types
      d_bond_nbond=enefunc%num_bond
      d_angle_nangle=enefunc%num_angle
      d_dihed_ndihe=enefunc%num_dihedral
      d_improper_nimp=enefunc%num_improper
      d_nb14_cutoff=enefunc%cutoffdist
      d_nb14_density=enefunc%table%density
      if(first) d_nb14_nb14_lj12=enefunc%nb14_lj12
      if(first) d_nb14_nb14_lj6=enefunc%nb14_lj6
      d_nb14_num_nb14_calc=enefunc%num_nb14_calc
      d_nb14_start_atom=domain%start_atom
      if(first) d_nb14_table_ene=enefunc%table%table_ene
      if(first) d_nb14_table_grad=enefunc%table%table_grad
      d_corr_cutoff=enefunc%cutoffdist
      d_corr_density=enefunc%table%density
      d_corr_num_nonb_excl=enefunc%num_nonb_excl
      d_corr_start_atom=domain%start_atom
      if(first) d_corr_table_decor=enefunc%table%table_decor
      if(first) d_corr_table_ecor=enefunc%table%table_ecor

    end if
    if(first.or.pairlist%univ_update/=0) then
      total=0
      do icell=1,cells
        h_offset_bond(icell)=total
        total=total+enefunc%num_bond(icell)
      end do
      if(total>capacity_bond) error stop 'GPU compact topology capacity changed'
      live_bond=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_bond(icell)
          slot=h_offset_bond(icell)+entry
          h_owner_bond(slot)=icell;h_local_bond(slot)=entry
          h_bond_bondlist(:,slot)=enefunc%bond_list(:,entry,icell)
          h_bond_fc(slot)=enefunc%bond_force_const(entry,icell)
          h_bond_r0(slot)=enefunc%bond_dist_min(entry,icell)
          h_bond_bond_pbc(slot)=enefunc%bond_pbc(entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_bond=h_offset_bond;owner_bond=h_owner_bond;local_bond=h_local_bond
      d_bond_bondlist=h_bond_bondlist
      d_bond_fc=h_bond_fc
      d_bond_r0=h_bond_r0
      d_bond_bond_pbc=h_bond_bond_pbc
      total=0
      do icell=1,cells
        h_offset_angle(icell)=total
        total=total+enefunc%num_angle(icell)
      end do
      if(total>capacity_angle) error stop 'GPU compact topology capacity changed'
      live_angle=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_angle(icell)
          slot=h_offset_angle(icell)+entry
          h_owner_angle(slot)=icell;h_local_angle(slot)=entry
          h_angle_anglelist(:,slot)=enefunc%angle_list(:,entry,icell)
          h_angle_fc(slot)=enefunc%angle_force_const(entry,icell)
          h_angle_theta0(slot)=enefunc%angle_theta_min(entry,icell)
          h_angle_fc_ub(slot)=enefunc%urey_force_const(entry,icell)
          h_angle_r0_ub(slot)=enefunc%urey_rmin(entry,icell)
          h_angle_angl_pbc(:,slot)=enefunc%angle_pbc(:,entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_angle=h_offset_angle;owner_angle=h_owner_angle;local_angle=h_local_angle
      d_angle_anglelist=h_angle_anglelist
      d_angle_fc=h_angle_fc
      d_angle_theta0=h_angle_theta0
      d_angle_fc_ub=h_angle_fc_ub
      d_angle_r0_ub=h_angle_r0_ub
      d_angle_angl_pbc=h_angle_angl_pbc
      total=0
      do icell=1,cells
        h_offset_dihed(icell)=total
        total=total+enefunc%num_dihedral(icell)
      end do
      if(total>capacity_dihed) error stop 'GPU compact topology capacity changed'
      live_dihed=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_dihedral(icell)
          slot=h_offset_dihed(icell)+entry
          h_owner_dihed(slot)=icell;h_local_dihed(slot)=entry
          h_dihed_dihelist(:,slot)=enefunc%dihe_list(:,entry,icell)
          h_dihed_fc(slot)=enefunc%dihe_force_const(entry,icell)
          h_dihed_nperiod(slot)=enefunc%dihe_periodicity(entry,icell)
          h_dihed_phase(slot)=enefunc%dihe_phase(entry,icell)
          h_dihed_dihe_pbc(:,slot)=enefunc%dihe_pbc(:,entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_dihed=h_offset_dihed;owner_dihed=h_owner_dihed;local_dihed=h_local_dihed
      d_dihed_dihelist=h_dihed_dihelist
      d_dihed_fc=h_dihed_fc
      d_dihed_nperiod=h_dihed_nperiod
      d_dihed_phase=h_dihed_phase
      d_dihed_dihe_pbc=h_dihed_dihe_pbc
      total=0
      do icell=1,cells
        h_offset_improper(icell)=total
        total=total+enefunc%num_improper(icell)
      end do
      if(total>capacity_improper) error stop 'GPU compact topology capacity changed'
      live_improper=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_improper(icell)
          slot=h_offset_improper(icell)+entry
          h_owner_improper(slot)=icell;h_local_improper(slot)=entry
          h_improper_imprlist(:,slot)=enefunc%impr_list(:,entry,icell)
          h_improper_fc(slot)=enefunc%impr_force_const(entry,icell)
          h_improper_phase(slot)=enefunc%impr_phase(entry,icell)
          h_improper_impr_pbc(:,slot)=enefunc%impr_pbc(:,entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_improper=h_offset_improper;owner_improper=h_owner_improper;local_improper=h_local_improper
      d_improper_imprlist=h_improper_imprlist
      d_improper_fc=h_improper_fc
      d_improper_phase=h_improper_phase
      d_improper_impr_pbc=h_improper_impr_pbc
      total=0
      do icell=1,cells
        h_offset_nb14(icell)=total
        total=total+enefunc%num_nb14_calc(icell)
      end do
      if(total>capacity_nb14) error stop 'GPU compact topology capacity changed'
      live_nb14=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_nb14_calc(icell)
          slot=h_offset_nb14(icell)+entry
          h_owner_nb14(slot)=icell;h_local_nb14(slot)=entry
          h_nb14_cell_move(:,slot)=domain%cell_move(:,enefunc%nb14_calc_list(2,entry,icell),enefunc%nb14_calc_list(1,entry,icell))
          h_nb14_nb14_calc_list(:,slot)=enefunc%nb14_calc_list(:,entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_nb14=h_offset_nb14;owner_nb14=h_owner_nb14;local_nb14=h_local_nb14
      d_nb14_cell_move=h_nb14_cell_move
      d_nb14_nb14_calc_list=h_nb14_nb14_calc_list
      total=0
      do icell=1,cells
        h_offset_corr(icell)=total
        total=total+enefunc%num_nonb_excl(icell)
      end do
      if(total>capacity_corr) error stop 'GPU compact topology capacity changed'
      live_corr=total
      !$omp parallel do private(entry,slot) schedule(static)
      do icell=1,cells
        do entry=1,enefunc%num_nonb_excl(icell)
          slot=h_offset_corr(icell)+entry
          h_owner_corr(slot)=icell;h_local_corr(slot)=entry
          h_corr_cell_move(:,slot)=domain%cell_move(:,enefunc%nonb_excl_list(2,entry,icell),enefunc%nonb_excl_list(1,entry,icell))
          h_corr_nonb_excl_list(:,slot)=enefunc%nonb_excl_list(:,entry,icell)
        end do
      end do
      !$omp end parallel do
      offset_corr=h_offset_corr;owner_corr=h_owner_corr;local_corr=h_local_corr
      d_corr_cell_move=h_corr_cell_move
      d_corr_nonb_excl_list=h_corr_nonb_excl_list
    end if
    snapshot_box=domain%system_size
    call checked(cudaMemcpyAsync(g_box,snapshot_box,3,cudaMemcpyHostToDevice,producer_stream))
    copied=step_export_force_coord(g_coord,producer_stream)
    if(.not.copied) then
      snapshot_coord=coord(:,1:width,1:cells)
      call checked(cudaMemcpyAsync(g_coord,snapshot_coord,size(g_coord),cudaMemcpyHostToDevice,producer_stream))
    end if
    if(coordinate_copy_c(g_pbc,atoms,4,producer_stream,1)==0) then
    snapshot_pbc=pbc(1:4*atoms,1:1,1:1)
    call checked(cudaMemcpyAsync(g_pbc,snapshot_pbc,size(g_pbc),cudaMemcpyHostToDevice,producer_stream))
    end if
    call checked(cudaMemsetAsync(g_force,0.0d0,size(g_force),producer_stream))
    call values_pack(before_values,vir,b,a,u,d,i,e,v)
    if(sum(enefunc%num_bond(1:cells))>0) call kernel_bond<<<(sum(enefunc%num_bond(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_bond<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())
    if(sum(enefunc%num_angle(1:cells))>0) call kernel_angle<<<(sum(enefunc%num_angle(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_angle<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())
    if(sum(enefunc%num_dihedral(1:cells))>0) call kernel_dihed<<<(sum(enefunc%num_dihedral(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_dihed<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())
    if(sum(enefunc%num_improper(1:cells))>0) call kernel_improper<<<(sum(enefunc%num_improper(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_improper<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())
    if(sum(enefunc%num_nb14_calc(1:cells))>0) call kernel_nb14<<<(sum(enefunc%num_nb14_calc(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_nb14<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())
    if(sum(enefunc%num_nonb_excl(1:cells))>0) call kernel_corr<<<(sum(enefunc%num_nonb_excl(1:cells))+127)/128,128,0,producer_stream>>>()
    call checked(cudaGetLastError())
    call reduce_corr<<<cells,64,0,producer_stream>>>()
    call checked(cudaGetLastError())

    call checked(cudaforSetDefaultStream(stream=previous_stream))
    first=.false.
  end subroutine

  subroutine bonded_finish(domain,vir,b,a,u,d,i,e,v)
    type(s_domain),intent(in)::domain
    real(dp),intent(inout)::vir(:,:,:),b(:),a(:),u(:),d(:),i(:),e(:),v(:)
    real(8)::sum_values(10,nthread),current(10,nthread),reference(10,nthread),ratio,maximum,atol
    integer::term,cell,worker,field
    integer(kind=cuda_stream_kind)::previous_stream
    previous_stream=cudaforGetDefaultStream()
    call checked(cudaforSetDefaultStream(stream=producer_stream))
    host_parts(:,:,1)=sum_bond
    host_parts(:,:,2)=sum_angle
    host_parts(:,:,3)=sum_dihed
    host_parts(:,:,4)=sum_improper
    host_parts(:,:,5)=sum_nb14
    host_parts(:,:,6)=sum_corr

    call checked(cudaStreamSynchronize(producer_stream))
    call checked(cudaforSetDefaultStream(stream=previous_stream))
    sum_values=0
    do term=1,6
      do cell=1,cells
        worker=mod(cell-1,nthread)+1
        sum_values(:,worker)=sum_values(:,worker)+host_parts(:,cell,term)
      end do
    end do
    if(.not.all(ieee_is_finite(sum_values))) error stop 'Nonfinite GPU bonded energy or virial'
    if(bonded_verify()) then
      call values_pack(current,vir,b,a,u,d,i,e,v)
      reference=current-before_values;maximum=0
      do field=1,10
        atol=1.0d-8
        if(field>=8) atol=1.0d-7
        ratio=maxval(abs(sum_values(field,:)-reference(field,:))/(atol+1.0d-10*abs(reference(field,:))))
        maximum=max(maximum,ratio)
      end do
      if(maximum>1) then
        print *, 'GENESIS_BONDED_VERIFY error ',maximum
        print *, 'CPU bonded totals ',sum(reference,dim=2)
        print *, 'GPU bonded totals ',sum(sum_values,dim=2)
        error stop 'GPU bonded energy or virial differs from CPU'
      end if
      write(*,'(a,es20.12)') 'GENESIS_BONDED_VERIFY: energy and virial pass ',maximum
    else
      b=b+sum_values(1,:);a=a+sum_values(2,:);u=u+sum_values(3,:);d=d+sum_values(4,:);i=i+sum_values(5,:)
      e=e+sum_values(6,:);v=v+sum_values(7,:)
      vir(1,1,:)=vir(1,1,:)+sum_values(8,:)
      vir(2,2,:)=vir(2,2,:)+sum_values(9,:)
      vir(3,3,:)=vir(3,3,:)+sum_values(10,:)
    end if
    call register_force(g_force)
  end subroutine

  attributes(global) subroutine kernel_bond()
    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc
    real(wp)                 :: ebond_temp, work(3), viri(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id, omp_get_thread_num
    integer                  :: list(3), pbc_int, k1, k2, k3

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_bond) return
    i=owner_bond(slot);ix=local_bond(slot)
    viri=0
    ebond_temp=0

        icel1 = g_id(1,d_bond_bondlist(1,slot))
        i1    = g_id(2,d_bond_bondlist(1,slot))
        icel2 = g_id(1,d_bond_bondlist(2,slot))
        i2    = g_id(2,d_bond_bondlist(2,slot))
        pbc_int = d_bond_bond_pbc(slot)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1

        ! bond energy: E=K[b-b0]^2
        !
        d12(1) = g_coord(1,i1,icel1) - g_coord(1,i2,icel2) &
               + g_box(1)*real(k1,wp)
        d12(2) = g_coord(2,i1,icel1) - g_coord(2,i2,icel2) &
               + g_box(2)*real(k2,wp)
        d12(3) = g_coord(3,i1,icel1) - g_coord(3,i2,icel2) &
               + g_box(3)*real(k3,wp)
        r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
        r_dif = r12 - d_bond_r0(slot)
        ebond_temp = ebond_temp + d_bond_fc(slot) * r_dif * r_dif

        ! gradient: dE/dX
        !
        cc_frc  = (2.0_wp * d_bond_fc(slot) * r_dif) / r12
        work(1) = cc_frc * d12(1)
        work(2) = cc_frc * d12(2)
        work(3) = cc_frc * d12(3)

        ! virial
        !
        viri(1) = viri(1) + d12(1)*work(1)  
        viri(2) = viri(2) + d12(2)*work(2)  
        viri(3) = viri(3) + d12(3)*work(3)  
  
        ! store force: F=-dE/dX
        !
    atomic_tmp=atomicadd(g_force(1,i1,icel1), - work(1))
    atomic_tmp=atomicadd(g_force(2,i1,icel1), - work(2))
    atomic_tmp=atomicadd(g_force(3,i1,icel1), - work(3))
    atomic_tmp=atomicadd(g_force(1,i2,icel2), + work(1))
    atomic_tmp=atomicadd(g_force(2,i2,icel2), + work(2))
    atomic_tmp=atomicadd(g_force(3,i2,icel2), + work(3))

    out_bond(:,slot)=0
    out_bond(1,slot)=ebond_temp
    out_bond(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_bond()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_bond(cell)+lane,offset_bond(cell)+d_bond_nbond(cell),64
      value=value+out_bond(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_bond(:,cell)=partial(:,1)
  end subroutine
  attributes(global) subroutine kernel_angle()
    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc_ub
    real(wp)                 :: eangle_temp, eurey_temp, work(9), viri(3)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: list(3), pbc_int, k1, k2, k3
    integer                  :: id, omp_get_thread_num

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_angle) return
    i=owner_angle(slot);ix=local_angle(slot)
    viri=0
    eangle_temp=0
    eurey_temp=0

        list(1) = d_angle_anglelist(1,slot)
        list(2) = d_angle_anglelist(2,slot)
        list(3) = d_angle_anglelist(3,slot)
        icel1   = g_id(1,list(1))
        i1      = g_id(2,list(1))
        icel2   = g_id(1,list(2))
        i2      = g_id(2,list(2))
        icel3   = g_id(1,list(3))
        i3      = g_id(2,list(3))

        ! angle energy: E=K[t-t0]^2
        !
        pbc_int = d_angle_angl_pbc(1,slot)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1
        d12(1) = g_coord(1,i1,icel1) - g_coord(1,i2,icel2) &
               + g_box(1)*real(k1,wp)
        d12(2) = g_coord(2,i1,icel1) - g_coord(2,i2,icel2) &
               + g_box(2)*real(k2,wp)
        d12(3) = g_coord(3,i1,icel1) - g_coord(3,i2,icel2) &
               + g_box(3)*real(k3,wp)

        pbc_int = d_angle_angl_pbc(2,slot)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1
        d32(1) = g_coord(1,i3,icel3) - g_coord(1,i2,icel2) &
               + g_box(1)*real(k1,wp)
        d32(2) = g_coord(2,i3,icel3) - g_coord(2,i2,icel2) &
               + g_box(2)*real(k2,wp)
        d32(3) = g_coord(3,i3,icel3) - g_coord(3,i2,icel2) &
               + g_box(3)*real(k3,wp)
        r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
        r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
        r12r32   = sqrt(r12_2*r32_2)

        inv_r12r32 = 1.0_wp / r12r32
        inv_r12_2  = 1.0_wp / r12_2
        inv_r32_2  = 1.0_wp / r32_2

        cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32
        cos_t  = min(1.0_wp, cos_t) 
        cos_t  = max(-1.0_wp, cos_t) 
        t123   = acos(cos_t)
        t_dif  = t123 - d_angle_theta0(slot)
        eangle_temp = eangle_temp + d_angle_fc(slot) * t_dif * t_dif
    
        ! gradient: dE/dX
        !
        sin_t  = sin(t123)
        sin_t  = max (EPS, sin_t)
        cc_frc = - (2.0_wp * d_angle_fc(slot) * t_dif) / sin_t
        cc_frc2 = cos_t * inv_r12_2
        cc_frc3 = cos_t * inv_r32_2
        work(1) = cc_frc * (d32(1) * inv_r12r32 - d12(1) * cc_frc2)
        work(2) = cc_frc * (d32(2) * inv_r12r32 - d12(2) * cc_frc2)
        work(3) = cc_frc * (d32(3) * inv_r12r32 - d12(3) * cc_frc2)
        work(4) = cc_frc * (d12(1) * inv_r12r32 - d32(1) * cc_frc3)
        work(5) = cc_frc * (d12(2) * inv_r12r32 - d32(2) * cc_frc3)
        work(6) = cc_frc * (d12(3) * inv_r12r32 - d32(3) * cc_frc3)

        viri(1) = viri(1) + d12(1)*work(1) + d32(1)*work(4)
        viri(2) = viri(2) + d12(2)*work(2) + d32(2)*work(5)
        viri(3) = viri(3) + d12(3)*work(3) + d32(3)*work(6)

        ! urey-bradley energy: E=K[r-r0]^2
        !
        if (abs(d_angle_fc_ub(slot)) > EPS) then

          pbc_int = d_angle_angl_pbc(3,slot)
          k3 = pbc_int / 9
          pbc_int = pbc_int - k3*9
          k2 = pbc_int / 3
          k1 = pbc_int - k2*3
          k1 = k1 - 1
          k2 = k2 - 1
          k3 = k3 - 1
          d13(1) = g_coord(1,i1,icel1) - g_coord(1,i3,icel3) &
                 + g_box(1)*real(k1,wp)
          d13(2) = g_coord(2,i1,icel1) - g_coord(2,i3,icel3) &
                 + g_box(2)*real(k2,wp)
          d13(3) = g_coord(3,i1,icel1) - g_coord(3,i3,icel3) &
                 + g_box(3)*real(k3,wp)
          r13    = sqrt(d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3))
          ub_dif = r13 - d_angle_r0_ub(slot)
          eurey_temp = eurey_temp + d_angle_fc_ub(slot) * ub_dif * ub_dif

          ! gradient: dE/dx
          !
          cc_frc_ub = (2.0_wp * d_angle_fc_ub(slot) * ub_dif) / r13
          work(7) = cc_frc_ub * d13(1)
          work(8) = cc_frc_ub * d13(2)
          work(9) = cc_frc_ub * d13(3)

          viri(1) = viri(1) + d13(1)*work(7)
          viri(2) = viri(2) + d13(2)*work(8)
          viri(3) = viri(3) + d13(3)*work(9)

        else

          work(7) = 0.0_wp
          work(8) = 0.0_wp
          work(9) = 0.0_wp

        end if

        ! store force
        !
    atomic_tmp=atomicadd(g_force(1,i1,icel1),-work(1)-work(7))
    atomic_tmp=atomicadd(g_force(2,i1,icel1),-work(2)-work(8))
    atomic_tmp=atomicadd(g_force(3,i1,icel1),-work(3)-work(9))
    atomic_tmp=atomicadd(g_force(1,i2,icel2),+work(1)+work(4))
    atomic_tmp=atomicadd(g_force(2,i2,icel2),+work(2)+work(5))
    atomic_tmp=atomicadd(g_force(3,i2,icel2),+work(3)+work(6))
    atomic_tmp=atomicadd(g_force(1,i3,icel3),-work(4)+work(7))
    atomic_tmp=atomicadd(g_force(2,i3,icel3),-work(5)+work(8))
    atomic_tmp=atomicadd(g_force(3,i3,icel3),-work(6)+work(9))

    out_angle(:,slot)=0
    out_angle(2,slot)=eangle_temp
    out_angle(3,slot)=eurey_temp
    out_angle(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_angle()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_angle(cell)+lane,offset_angle(cell)+d_angle_nangle(cell),64
      value=value+out_angle(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_angle(:,cell)=partial(:,1)
  end subroutine
  attributes(global) subroutine kernel_dihed()
    ! local variables
    integer                  :: i, j, k, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, edihe_temp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(9), v(3), cwork(3,4)
    real(wp)                 :: work(9), viri(3)
    integer                  :: omp_get_thread_num, list(4), pbc_index(3)

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_dihed) return
    i=owner_dihed(slot);ix=local_dihed(slot)
    viri=0
    edihe_temp=0

        list(1) = d_dihed_dihelist(1,slot)
        list(2) = d_dihed_dihelist(2,slot)
        list(3) = d_dihed_dihelist(3,slot)
        list(4) = d_dihed_dihelist(4,slot)

        icel1 = g_id(1,list(1))
        i1    = g_id(2,list(1))
        icel2 = g_id(1,list(2))
        i2    = g_id(2,list(2))
        icel3 = g_id(1,list(3))
        i3    = g_id(2,list(3))
        icel4 = g_id(1,list(4))
        i4    = g_id(2,list(4))

        cwork(1,1) = g_coord(1,i1,icel1)
        cwork(2,1) = g_coord(2,i1,icel1)
        cwork(3,1) = g_coord(3,i1,icel1)
        cwork(1,2) = g_coord(1,i2,icel2)
        cwork(2,2) = g_coord(2,i2,icel2)
        cwork(3,2) = g_coord(3,i2,icel2)
        cwork(1,3) = g_coord(1,i3,icel3)
        cwork(2,3) = g_coord(2,i3,icel3)
        cwork(3,3) = g_coord(3,i3,icel3)
        cwork(1,4) = g_coord(1,i4,icel4)
        cwork(2,4) = g_coord(2,i4,icel4)
        cwork(3,4) = g_coord(3,i4,icel4)

        pbc_index(1) = d_dihed_dihe_pbc(1,slot)
        pbc_index(2) = d_dihed_dihe_pbc(2,slot)
        pbc_index(3) = d_dihed_dihe_pbc(3,slot)
        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral_2( &
                              aindex, pbc_index, cwork, g_box,  &
                              cos_dih, sin_dih, grad, v)
        cosnt = 1.0_wp
        sinnt = 0.0_wp
        krot = 0
        nrot = d_dihed_nperiod(slot)
        if (g_notation > 0) &
        nrot = mod(d_dihed_nperiod(slot), g_notation)
        do while (krot < nrot)
          tmp   = cosnt*cos_dih - sinnt*sin_dih
          sinnt = sinnt*cos_dih + cosnt*sin_dih
          cosnt = tmp
          krot = krot+1
        end do

        cospha = cos(d_dihed_phase(slot))
        sinpha = sin(d_dihed_phase(slot))
        edihe_temp = edihe_temp + d_dihed_fc(slot)   &
                                 *(1.0_wp+cospha*cosnt+sinnt*sinpha)

        grad_coef = d_dihed_fc(slot) * real(nrot,wp)    &
                              * (cospha*sinnt-cosnt*sinpha)
        work(1) = grad_coef*grad(1)
        work(2) = grad_coef*grad(2)
        work(3) = grad_coef*grad(3)
        work(4) = grad_coef*grad(4)
        work(5) = grad_coef*grad(5)
        work(6) = grad_coef*grad(6)
        work(7) = grad_coef*grad(7)
        work(8) = grad_coef*grad(8)
        work(9) = grad_coef*grad(9)

        viri(1) = viri(1) + grad_coef*v(1)
        viri(2) = viri(2) + grad_coef*v(2)
        viri(3) = viri(3) + grad_coef*v(3)

    atomic_tmp=atomicadd(g_force(1,i1,icel1), - work(1))
    atomic_tmp=atomicadd(g_force(2,i1,icel1), - work(2))
    atomic_tmp=atomicadd(g_force(3,i1,icel1), - work(3))
    atomic_tmp=atomicadd(g_force(1,i2,icel2),                                  + work(1) - work(4))
    atomic_tmp=atomicadd(g_force(2,i2,icel2),                                  + work(2) - work(5))
    atomic_tmp=atomicadd(g_force(3,i2,icel2),                                  + work(3) - work(6))
    atomic_tmp=atomicadd(g_force(1,i3,icel3),                                  + work(4) + work(7))
    atomic_tmp=atomicadd(g_force(2,i3,icel3),                                  + work(5) + work(8))
    atomic_tmp=atomicadd(g_force(3,i3,icel3),                                  + work(6) + work(9))
    atomic_tmp=atomicadd(g_force(1,i4,icel4), - work(7))
    atomic_tmp=atomicadd(g_force(2,i4,icel4), - work(8))
    atomic_tmp=atomicadd(g_force(3,i4,icel4), - work(9))

    out_dihed(:,slot)=0
    out_dihed(4,slot)=edihe_temp
    out_dihed(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_dihed()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_dihed(cell)+lane,offset_dihed(cell)+d_dihed_ndihe(cell),64
      value=value+out_dihed(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_dihed(:,cell)=partial(:,1)
  end subroutine
  attributes(global) subroutine kernel_improper()
    ! local variables
    integer                  :: i, ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(4), pbc_index(3)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, eimp_temp
    real(wp)                 :: cosdif, sindif, diffphi
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(9), v(3), cwork(3,4)
    real(wp)                 :: work(9), viri(3)
    integer                  :: id, omp_get_thread_num

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_improper) return
    i=owner_improper(slot);ix=local_improper(slot)
    viri=0
    eimp_temp=0

        icel1 = g_id(1,d_improper_imprlist(1,slot))
        i1    = g_id(2,d_improper_imprlist(1,slot))
        icel2 = g_id(1,d_improper_imprlist(2,slot))
        i2    = g_id(2,d_improper_imprlist(2,slot))
        icel3 = g_id(1,d_improper_imprlist(3,slot))
        i3    = g_id(2,d_improper_imprlist(3,slot))
        icel4 = g_id(1,d_improper_imprlist(4,slot))
        i4    = g_id(2,d_improper_imprlist(4,slot))

        cwork(1,1) = g_coord(1,i1,icel1)
        cwork(2,1) = g_coord(2,i1,icel1)
        cwork(3,1) = g_coord(3,i1,icel1)
        cwork(1,2) = g_coord(1,i2,icel2)
        cwork(2,2) = g_coord(2,i2,icel2)
        cwork(3,2) = g_coord(3,i2,icel2)
        cwork(1,3) = g_coord(1,i3,icel3)
        cwork(2,3) = g_coord(2,i3,icel3)
        cwork(3,3) = g_coord(3,i3,icel3)
        cwork(1,4) = g_coord(1,i4,icel4)
        cwork(2,4) = g_coord(2,i4,icel4)
        cwork(3,4) = g_coord(3,i4,icel4)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        pbc_index(1) = d_improper_impr_pbc(1,slot)
        pbc_index(2) = d_improper_impr_pbc(2,slot)
        pbc_index(3) = d_improper_impr_pbc(3,slot)

        call calculate_dihedral_2( &
                              aindex, pbc_index, cwork, g_box, &
                              cos_dih, sin_dih, grad, v)

        cospha = cos(d_improper_phase(slot))
        sinpha = sin(d_improper_phase(slot))

        cosdif = cos_dih*cospha + sin_dih*sinpha
        sindif = cos_dih*sinpha - sin_dih*cospha

        if (cosdif > 1.0E-1_wp) then
          diffphi = asin(sindif)
        else
          diffphi = sign(1.0_wp,sindif)*acos(cosdif)
        end if
        eimp_temp = eimp_temp + d_improper_fc(slot)*diffphi*diffphi
        grad_coef = 2.0_wp*d_improper_fc(slot)*diffphi

        work(1) = grad_coef*grad(1)
        work(2) = grad_coef*grad(2)
        work(3) = grad_coef*grad(3)
        work(4) = grad_coef*grad(4)
        work(5) = grad_coef*grad(5)
        work(6) = grad_coef*grad(6)
        work(7) = grad_coef*grad(7)
        work(8) = grad_coef*grad(8)
        work(9) = grad_coef*grad(9)

        viri(1) = viri(1) + grad_coef*v(1)
        viri(2) = viri(2) + grad_coef*v(2)
        viri(3) = viri(3) + grad_coef*v(3)

    atomic_tmp=atomicadd(g_force(1,i1,icel1), - work(1))
    atomic_tmp=atomicadd(g_force(2,i1,icel1), - work(2))
    atomic_tmp=atomicadd(g_force(3,i1,icel1), - work(3))
    atomic_tmp=atomicadd(g_force(1,i2,icel2), + work(1) - work(4))
    atomic_tmp=atomicadd(g_force(2,i2,icel2), + work(2) - work(5))
    atomic_tmp=atomicadd(g_force(3,i2,icel2), + work(3) - work(6))
    atomic_tmp=atomicadd(g_force(1,i3,icel3), + work(4) + work(7))
    atomic_tmp=atomicadd(g_force(2,i3,icel3), + work(5) + work(8))
    atomic_tmp=atomicadd(g_force(3,i3,icel3), + work(6) + work(9))
    atomic_tmp=atomicadd(g_force(1,i4,icel4), - work(7))
    atomic_tmp=atomicadd(g_force(2,i4,icel4), - work(8))
    atomic_tmp=atomicadd(g_force(3,i4,icel4), - work(9))

    out_improper(:,slot)=0
    out_improper(5,slot)=eimp_temp
    out_improper(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_improper()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_improper(cell)+lane,offset_improper(cell)+d_improper_nimp(cell),64
      value=value+out_improper(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_improper(:,cell)=partial(:,1)
  end subroutine
  attributes(global) subroutine kernel_nb14()
    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: R, table(6), iqtmp, jqtmp
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(3), viri(3)
    real(wp)                 :: evdw_temp, elec_temp
    integer                  :: i, ix, iy, j, k, ij, icel, jcel, m, L, L1
    integer                  :: iatmcls, jatmcls, num_atom
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_nb14) return
    ij=owner_nb14(slot);k=local_nb14(slot)
    viri=0
    elec_temp=0
    evdw_temp=0
    cutoff2=d_nb14_cutoff*d_nb14_cutoff

        icel  = d_nb14_nb14_calc_list(1,slot) 
        jcel  = d_nb14_nb14_calc_list(2,slot) 
        ix    = d_nb14_nb14_calc_list(3,slot)
        iy    = d_nb14_nb14_calc_list(4,slot)
        i     = d_nb14_start_atom(icel) + ix
        j     = d_nb14_start_atom(jcel) + iy
        iatmcls = g_class(i)
        jatmcls = g_class(j)

        ! compute distance
        !
        dij(1) = g_pbc(           i,1,1) - g_pbc(           j,1,1)
        dij(2) = g_pbc(  g_atoms+i,1,1) - g_pbc(  g_atoms+j,1,1)
        dij(3) = g_pbc(2*g_atoms+i,1,1) - g_pbc(2*g_atoms+j,1,1)
        dij(1) = dij(1) + d_nb14_cell_move(1,slot)*g_box(1)
        dij(2) = dij(2) + d_nb14_cell_move(2,slot)*g_box(2)
        dij(3) = dij(3) + d_nb14_cell_move(3,slot)*g_box(3)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        iqtmp  = g_pbc(3*g_atoms+i,1,1)
        jqtmp  = g_pbc(3*g_atoms+j,1,1) 

        ! energy and gradient
        !
        rij2  = cutoff2*d_nb14_density/rij2
        L     = int(rij2)
        R     = rij2 - L
        L1    = 3*L - 2

        lj6   = d_nb14_nb14_lj6 (iatmcls,jatmcls)
        lj12  = d_nb14_nb14_lj12(iatmcls,jatmcls)

        term_lj12 = d_nb14_table_ene(L1)   + R*(d_nb14_table_ene(L1+3)-d_nb14_table_ene(L1)  )
        term_lj6  = d_nb14_table_ene(L1+1) + R*(d_nb14_table_ene(L1+4)-d_nb14_table_ene(L1+1))
        term_elec = d_nb14_table_ene(L1+2) + R*(d_nb14_table_ene(L1+5)-d_nb14_table_ene(L1+2))
        evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
        elec_temp = elec_temp + iqtmp*jqtmp*term_elec

        term_lj12 = d_nb14_table_grad(L1)   + R*(d_nb14_table_grad(L1+3)-d_nb14_table_grad(L1)  )
        term_lj6  = d_nb14_table_grad(L1+1) + R*(d_nb14_table_grad(L1+4)-d_nb14_table_grad(L1+1))
        term_elec = d_nb14_table_grad(L1+2) + R*(d_nb14_table_grad(L1+5)-d_nb14_table_grad(L1+2))
        grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec

        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

    atomic_tmp=atomicadd(g_force(1,ix,icel), - work(1))
    atomic_tmp=atomicadd(g_force(2,ix,icel), - work(2))
    atomic_tmp=atomicadd(g_force(3,ix,icel), - work(3))
    atomic_tmp=atomicadd(g_force(1,iy,jcel), + work(1))
    atomic_tmp=atomicadd(g_force(2,iy,jcel), + work(2))
    atomic_tmp=atomicadd(g_force(3,iy,jcel), + work(3))

    out_nb14(:,slot)=0
    out_nb14(6,slot)=elec_temp
    out_nb14(7,slot)=evdw_temp
    out_nb14(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_nb14()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_nb14(cell)+lane,offset_nb14(cell)+d_nb14_num_nb14_calc(cell),64
      value=value+out_nb14(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_nb14(:,cell)=partial(:,1)
  end subroutine
  attributes(global) subroutine kernel_corr()
    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), viri(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef, elec_temp
    real(wp)                 :: iqtmp, jqtmp
    integer                  :: i, ix, iy, j, k, ij, L, icel, jcel
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local, num_atom

    integer::slot
    real(8)::atomic_tmp
    slot=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(slot>live_corr) return
    ij=owner_corr(slot);k=local_corr(slot)
    viri=0
    elec_temp=0
    cutoff2=d_corr_cutoff*d_corr_cutoff

        icel  = d_corr_nonb_excl_list(1,slot)
        jcel  = d_corr_nonb_excl_list(2,slot)
        ix    = d_corr_nonb_excl_list(3,slot)
        iy    = d_corr_nonb_excl_list(4,slot)
        i     = d_corr_start_atom(icel) + ix
        j     = d_corr_start_atom(jcel) + iy

        ! compute distance
        !
        dij(1) = g_pbc(           i,1,1) - g_pbc(           j,1,1)
        dij(2) = g_pbc(  g_atoms+i,1,1) - g_pbc(  g_atoms+j,1,1)
        dij(3) = g_pbc(2*g_atoms+i,1,1) - g_pbc(2*g_atoms+j,1,1)
        dij(1) = dij(1) + d_corr_cell_move(1,slot)*g_box(1)
        dij(2) = dij(2) + d_corr_cell_move(2,slot)*g_box(2)
        dij(3) = dij(3) + d_corr_cell_move(3,slot)*g_box(3)

        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        iqtmp = g_pbc(3*g_atoms+i,1,1)
        jqtmp = g_pbc(3*g_atoms+j,1,1)

        ! energy and gradient
        !
        rij2 = cutoff2*d_corr_density/rij2
        L    = int(rij2)
        R    = rij2 - L

        term_elec = d_corr_table_ecor(L) + R*(d_corr_table_ecor(L+1)-d_corr_table_ecor(L))
        ccr       = iqtmp * jqtmp * term_elec
        elec_temp = elec_temp + ccr

        term_elec = d_corr_table_decor(L) + R*(d_corr_table_decor(L+1)-d_corr_table_decor(L))
        coef      = iqtmp * jqtmp * term_elec

        work(1) = coef*dij(1)
        work(2) = coef*dij(2)
        work(3) = coef*dij(3)

        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

    atomic_tmp=atomicadd(g_force(1,ix,icel), - work(1))
    atomic_tmp=atomicadd(g_force(2,ix,icel), - work(2))
    atomic_tmp=atomicadd(g_force(3,ix,icel), - work(3))
    atomic_tmp=atomicadd(g_force(1,iy,jcel), + work(1))
    atomic_tmp=atomicadd(g_force(2,iy,jcel), + work(2))
    atomic_tmp=atomicadd(g_force(3,iy,jcel), + work(3))

    out_corr(:,slot)=0
    out_corr(6,slot)=elec_temp
    out_corr(8:10,slot)=-viri
  end subroutine
  attributes(global) subroutine reduce_corr()
    real(8),shared::partial(10,64)
    real(8)::value(10)
    integer::cell,item,k,stride,lane
    cell=blockIdx%x;lane=threadIdx%x;value=0
    do item=offset_corr(cell)+lane,offset_corr(cell)+d_corr_num_nonb_excl(cell),64
      value=value+out_corr(:,item)
    end do
    partial(:,lane)=value
    call syncthreads()
    stride=32
    do while(stride>0)
      if(lane<=stride) partial(:,lane)=partial(:,lane)+partial(:,lane+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(lane==1) sum_corr(:,cell)=partial(:,1)
  end subroutine
  attributes(device) subroutine calculate_dihedral_2(aindex, pbc_index, coord, box_size, &
                                  cos_dih, sin_dih, grad, v)

    ! formal arguments
    integer,   intent(in)    :: aindex(4)
    integer,   intent(in)    :: pbc_index(3)
    real(wp),  intent(in)    :: coord(3,4)
    real(wp),  intent(in)    :: box_size(3)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(9)
    real(wp),  intent(inout) :: v(3)

    ! local variable
    integer                  :: pbc_int, k1, k2, k3
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl


    pbc_int = pbc_index(1)
    k3 = pbc_int / 9
    pbc_int = pbc_int - k3*9
    k2 = pbc_int / 3
    k1 = pbc_int - k2*3
    k1 = k1 - 1
    k2 = k2 - 1
    k3 = k3 - 1
    dij(1) = coord(1,aindex(1)) - coord(1,aindex(2)) &
           + box_size(1)*real(k1,wp)
    dij(2) = coord(2,aindex(1)) - coord(2,aindex(2)) &
           + box_size(2)*real(k2,wp)
    dij(3) = coord(3,aindex(1)) - coord(3,aindex(2)) &
           + box_size(3)*real(k3,wp)

    pbc_int = pbc_index(2)
    k3 = pbc_int / 9
    pbc_int = pbc_int - k3*9
    k2 = pbc_int / 3
    k1 = pbc_int - k2*3
    k1 = k1 - 1
    k2 = k2 - 1
    k3 = k3 - 1
    djk(1) = coord(1,aindex(2)) - coord(1,aindex(3)) &
           + box_size(1)*real(k1,wp)
    djk(2) = coord(2,aindex(2)) - coord(2,aindex(3)) &
           + box_size(2)*real(k2,wp)
    djk(3) = coord(3,aindex(2)) - coord(3,aindex(3)) &
           + box_size(3)*real(k3,wp)

    pbc_int = pbc_index(3)
    k3 = pbc_int / 9
    pbc_int = pbc_int - k3*9
    k2 = pbc_int / 3
    k1 = pbc_int - k2*3
    k1 = k1 - 1
    k2 = k2 - 1
    k3 = k3 - 1
    dlk(1) = coord(1,aindex(4)) - coord(1,aindex(3)) &
           + box_size(1)*real(k1,wp)
    dlk(2) = coord(2,aindex(4)) - coord(2,aindex(3)) &
           + box_size(2)*real(k2,wp)
    dlk(3) = coord(3,aindex(4)) - coord(3,aindex(3)) &
           + box_size(3)*real(k3,wp)

    aijk(1) = dij(2)*djk(3) - dij(3)*djk(2)
    aijk(2) = dij(3)*djk(1) - dij(1)*djk(3)
    aijk(3) = dij(1)*djk(2) - dij(2)*djk(1)

    ajkl(1) = dlk(2)*djk(3) - dlk(3)*djk(2)
    ajkl(2) = dlk(3)*djk(1) - dlk(1)*djk(3)
    ajkl(3) = dlk(1)*djk(2) - dlk(2)*djk(1)

    raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
    rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)

    inv_raijk2  = 1.0_wp / raijk2
    inv_rajkl2  = 1.0_wp / rajkl2

    inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

    cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
    rjk     = sqrt(djk(1)*djk(1) + djk(2)*djk(2) + djk(3)*djk(3))
    inv_rjk = 1.0_wp/rjk

    tmp(1)  = aijk(1)*dlk(1) + aijk(2)*dlk(2) + aijk(3)*dlk(3)
    sin_dih = tmp(1) * rjk * inv_raijkl

    dotpro_ijk = dij(1)*djk(1) + dij(2)*djk(2) + dij(3)*djk(3)
    dotpro_jkl = djk(1)*dlk(1) + djk(2)*dlk(2) + djk(3)*dlk(3)

    tmp(1) = rjk*inv_raijk2
    tmp(2) = rjk*inv_rajkl2

    tmp(3) =  dotpro_ijk*inv_raijk2*inv_rjk
    tmp(4) =  dotpro_jkl*inv_rajkl2*inv_rjk

    grad(1) =  tmp(1)*aijk(1)
    grad(2) =  tmp(1)*aijk(2)
    grad(3) =  tmp(1)*aijk(3)
    grad(4) = -tmp(3)*aijk(1) + tmp(4)*ajkl(1)
    grad(5) = -tmp(3)*aijk(2) + tmp(4)*ajkl(2)
    grad(6) = -tmp(3)*aijk(3) + tmp(4)*ajkl(3)
    grad(7) =                 - tmp(2)*ajkl(1)
    grad(8) =                 - tmp(2)*ajkl(2)
    grad(9) =                 - tmp(2)*ajkl(3)

    v(1) = grad(1)*dij(1) + grad(4)*djk(1) + grad(7)*dlk(1)
    v(2) = grad(2)*dij(2) + grad(5)*djk(2) + grad(8)*dlk(2)
    v(3) = grad(3)*dij(3) + grad(6)*djk(3) + grad(9)*dlk(3)

    return

  end subroutine calculate_dihedral_2
end module
