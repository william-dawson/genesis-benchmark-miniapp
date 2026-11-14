#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module genesis_step_state_mod
  use iso_c_binding,only:c_int,c_double
  use cudafor
  use force_assembly_mod,only:force_resident_requested
  use sp_assign_velocity_mod,only:stop_translation
  use genesis_thermo_mod
  use coordinate_pack_mod
  use water_migration_mod,only:water_latch_c
  use hgr_migration_mod,only:hgr_latch_c
  use constants_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use sp_dynamics_str_mod
  use sp_ensemble_str_mod
  use sp_dynvars_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod, only: TIP3
  use mpi_parallel_mod
  use, intrinsic :: ieee_arithmetic, only:ieee_is_finite,ieee_value,ieee_quiet_nan
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
  implicit none
  private
  public :: step_thermo_candidate,step_thermo_begin,step_thermo_scale,step_thermo_stats
  public :: step_export_force_coord
  public :: step_stop_translation
  public :: step_enabled, step_verify, step_candidate2, step_run2, step_sync, step_forget
  public :: step_check_begin, step_check_reference, step_check_finish
  public :: step_candidate1,step_run1,step_check1_finish,step_coords,step_output_due,step_force2
  logical,save::thermo_pending=.false.
  logical,save::work_ready=.true.
  logical, save :: resident=.false.
  logical, save :: maps_current=.false.
  integer, save :: cells=0, ld=0, workers=0
  real(8), device, allocatable :: dx(:,:,:),dv(:,:,:),dxref(:,:,:),dvref(:,:,:),dwork(:,:,:)
  real(8), device, allocatable :: dfs(:,:,:),dfl(:,:,:),dmass(:,:),dimass(:,:),dcm(:,:)
  integer, device, allocatable :: dnatom(:),dnwater(:),dwater(:,:,:),dhlocal(:,:),dhlist(:,:,:,:)
  real(8), device, allocatable :: dhdist(:,:,:,:)
  integer, device :: failed
  integer,device,allocatable :: dnsolute(:)
  real(8),device,allocatable :: dwater_vir(:,:,:,:),dsolute_vir(:,:,:,:),dvir(:,:,:)
  real(8),allocatable,pinned :: host_vir(:,:,:)
  real(8),constant :: dscale(3)
  integer, constant :: dcells,diter,dconnect
  real(8), constant :: dtolerance,dwater_rOH,dwater_rHH,dwater_massO,dwater_massH
  real(8), constant :: dcenter(8),dbmoment(3)
  real(8), allocatable, pinned :: host_state(:,:,:,:),host_cm(:,:)
  real(8), allocatable :: saved(:,:,:,:),expected(:,:,:,:)
  real(8) :: saved_box(3),expected_box(3)
  type(s_boundary) :: saved_boundary,expected_boundary
  real(8),allocatable :: saved_trans(:,:,:)
  real(8) :: saved_system(3)
  interface
    subroutine poison_work_c(first,second,width,count) bind(C,name='genesis_poison_work')
      import c_double,c_int
      real(c_double),device::first(*),second(*)
      integer(c_int),value::width,count
    end subroutine
    subroutine remove_translation_c(velocity,counts,width,count,x,y,z) bind(C,name='genesis_stop_translation')
      import c_double,c_int
      real(c_double),device::velocity(*)
      integer(c_int),device,intent(in)::counts(*)
      integer(c_int),value::width,count
      real(c_double),value::x,y,z
    end subroutine
    integer(c_int) function assembly_copy(values,width,count) bind(C,name='genesis_force_copy')
      use iso_c_binding,only:c_int,c_double
      real(c_double),device::values(*)
      integer(c_int),value::width,count
    end function
  end interface
contains
  logical function step_thermo_candidate(dynamics,istep,istart,ensemble,domain,constraints)
    type(s_dynamics),intent(in)::dynamics
    integer,intent(in)::istep,istart
    type(s_ensemble),intent(in)::ensemble
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    step_thermo_candidate=.false.
    if(.not.thermo_enabled()) return
    if(.not.step_candidate2(ensemble,domain,constraints)) return
    if(istep<=istart.or.dynamics%thermo_period<=0) return
    if(istep>dynamics%nsteps) return
    if(dynamics%thermo_period/=dynamics%baro_period) return
    step_thermo_candidate=mod(istep-1,dynamics%thermo_period)==0
  end function

  subroutine step_thermo_begin(domain,constraints)
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    call prepare(domain,constraints)
    call save_reference<<<(ld*cells+255)/256,256>>>(ld,cells)
    call checked(cudaGetLastError())
  end subroutine

  subroutine step_thermo_scale(domain,constraints,factor)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(8),intent(in)::factor(3)
    call thermo_scale_reference(domain,constraints,dv,factor)
    dscale=factor
    call scale_groups<<<dim3(cells,(ld+63)/64,constraints%connect+2),64>>>(1)
    call checked(cudaGetLastError())
    call thermo_scale_check(domain,dv)
  end subroutine

  subroutine step_thermo_stats(domain,constraints,half,values)
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    real(8),intent(in)::half
    real(8),intent(out)::values(9)
    integer::cell,n
    call thermo_statistics(domain,constraints,dv,dxref,dfl,dfs,dmass,dimass, &
      dnatom,dnsolute,dnwater,dwater,dhlocal,dhlist,half,values)
    thermo_pending=.true.
    if(thermo_verify()) then
      do cell=1,cells
        n=domain%num_atom(cell)
        domain%velocity(:,1:n,cell)=ieee_value(0.0_wip,ieee_quiet_nan)
        domain%force(:,1:n,cell)=ieee_value(0.0_wip,ieee_quiet_nan)
        domain%velocity_half(:,1:n,cell)=ieee_value(0.0_wip,ieee_quiet_nan)
      end do
      write(*,'(a)') 'GENESIS_THERMO_VERIFY: CPU mirrors poisoned'
    end if
  end subroutine

  logical function state_work_enabled()
    integer,save::selected=-1
    integer::status,code
    character(16)::value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_STATE_WORK',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0.or.selected<0.or.selected>1) error stop 'GENESIS_STATE_WORK must be 0 or 1'
      end if
    end if
    state_work_enabled=selected==1
  end function

  logical function state_verify_enabled()
    integer::status
    character(16)::value
    call get_environment_variable('GENESIS_STATE_VERIFY',value,status=status)
    state_verify_enabled=status==0.and.trim(value)=='1'
  end function

  logical function step_stop_translation(domain,stop_trans) result(handled)
    type(s_domain),intent(inout)::domain
    logical,intent(in)::stop_trans
    integer,save::selected=-1
    integer::status,code,cell,atom
    real(8)::cm(4),ratio,maximum
    real(8),allocatable::reference(:,:,:)
    logical::verify
    character(16)::value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_STATE_GPU',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0.or.selected<0.or.selected>1) error stop 'GENESIS_STATE_GPU must be 0 or 1'
      end if
    end if
    handled=selected==1.and.resident.and.nproc_country==1
    if(.not.handled.or..not.stop_trans) return
    call get_environment_variable('GENESIS_STATE_VERIFY',value,status=status)
    verify=status==0.and.trim(value)=='1'
    if(verify) then
      call step_sync(domain)
      allocate(reference,mold=domain%velocity)
      reference=domain%velocity
      call stop_translation(domain%num_cell_local,domain%num_atom,stop_trans,domain%mass,reference)
    end if
    call center_momentum<<<cells,128>>>(cells)
    call checked(cudaGetLastError())
    host_cm=dcm;cm=0
    do cell=1,cells
      cm=cm+host_cm(1:4,cell)
    end do
    if(cm(4)<=0) error stop 'GPU translation removal has invalid total mass'
    cm(1:3)=cm(1:3)/cm(4)
    call remove_translation_c(dv,dnatom,ld,cells,cm(1),cm(2),cm(3))
    if(verify) then
      host_state(:,:,:,1)=dv;maximum=0
      do cell=1,cells
        do atom=1,domain%num_atom(cell)
          if(.not.all(ieee_is_finite(host_state(:,atom,cell,1)))) error stop 'Nonfinite GPU translation velocity'
          ratio=maxval(abs(host_state(:,atom,cell,1)-reference(:,atom,cell))/ &
              (1.0d-10+1.0d-12*abs(reference(:,atom,cell))))
          maximum=max(maximum,ratio)
        end do
      end do
      if(maximum>1) error stop 'GPU translation removal differs from CPU'
      domain%velocity=ieee_value(0.0_wip,ieee_quiet_nan)
      write(*,'(a,es20.12)') 'GENESIS_STATE_VERIFY: translation pass ',maximum
    end if
  end function

  logical function step_export_force_coord(values,stream)
    real(8),device,intent(out)::values(:,:,:)
    integer(kind=cuda_stream_kind),intent(in)::stream
    step_export_force_coord=resident
    if(resident) call checked(cudaMemcpyAsync(values,dx,size(dx),cudaMemcpyDeviceToDevice,stream))
  end function
  logical function step_enabled()
    integer, save :: selected=-1
    integer :: status,code
    character(16) :: value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_STEP_GPU',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0) error stop 'Invalid GENESIS_STEP_GPU'
        if(selected<0.or.selected>1) error stop 'GENESIS_STEP_GPU must be 0 or 1'
      end if
    end if
    step_enabled=selected==1
  end function

  logical function step_output_due(dynamics,istep)
    type(s_dynamics),intent(in)::dynamics
    integer,intent(in)::istep
    step_output_due=.false.
    if(istep<=0) return
    if(dynamics%crdout_period>0) step_output_due=mod(istep,dynamics%crdout_period)==0
    if(dynamics%velout_period>0) step_output_due=step_output_due.or.mod(istep,dynamics%velout_period)==0
    if(dynamics%rstout_period>0) step_output_due=step_output_due.or.mod(istep,dynamics%rstout_period)==0
  end function

  logical function step_verify()
    character(16) :: value
    integer :: status
    call get_environment_variable('GENESIS_STEP_VERIFY',value,status=status)
    step_verify=status==0.and.(trim(value)=='1'.or.trim(value)=='2')
  end function

  logical function step_diagnose()
    character(16)::value
    integer::status
    call get_environment_variable('GENESIS_STEP_VERIFY',value,status=status)
    step_diagnose=status==0.and.trim(value)=='2'
  end function

  logical function step_candidate2(ensemble,domain,constraints)
    type(s_ensemble), intent(in) :: ensemble
    type(s_domain), intent(in) :: domain
    type(s_constraints), intent(in) :: constraints
    step_candidate2=step_enabled().and.ensemble%ensemble==EnsembleNPT.and.ensemble%group_tp &
        .and.constraints%water_type==TIP3.and.constraints%rigid_bond &
        .and.constraints%connect<=7.and.domain%num_cell_boundary==0.and.nproc_country==1 &
        .and.ensemble%tpcontrol==TpcontrolBussi
#if defined(_SINGLE) || defined(_MIXED)
    step_candidate2=.false.
#endif
  end function

  logical function step_candidate1(dynamics,istep,ensemble,domain,constraints)
    type(s_dynamics),intent(in)::dynamics
    integer,intent(in)::istep
    type(s_ensemble),intent(in)::ensemble
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    step_candidate1=step_candidate2(ensemble,domain,constraints).and.istep>1
    step_candidate1=step_candidate1.and.mod(istep-1,dynamics%thermo_period)/=0 &
                    .and.mod(istep-1,dynamics%baro_period)/=0
  end function

  subroutine checked(status)
    integer,intent(in)::status
    if(status/=cudaSuccess) then
      print *, 'GENESIS_STEP CUDA failure: ',status
      error stop 'GPU timestep runtime error'
    end if
  end subroutine

  subroutine step_forget(remap)
    logical,optional,intent(in)::remap
    resident=.false.
    call coordinate_invalidate_c()
    if(present(remap)) then
      if(remap) maps_current=.false.
    end if
  end subroutine

  subroutine pack_state(domain,copy_short,copy_work)
    type(s_domain), intent(in) :: domain
    logical,intent(in)::copy_short,copy_work
    integer :: cell,atom
    !$omp parallel do private(atom) schedule(static)
    do cell=1,cells
      do atom=1,domain%num_atom(cell)
        host_state(:,atom,cell,1)=domain%velocity(:,atom,cell)
        host_state(:,atom,cell,2)=domain%coord(:,atom,cell)
        if(copy_work) host_state(:,atom,cell,3)=domain%coord_ref(:,atom,cell)
        if(copy_work) host_state(:,atom,cell,4)=domain%velocity_full(:,atom,cell)
        host_state(:,atom,cell,5)=domain%force_long(:,atom,cell)
        if(copy_short) host_state(:,atom,cell,6)=domain%force_short(:,atom,cell)
        host_state(:,atom,cell,7)=domain%velocity_ref(:,atom,cell)
      end do
    end do
    !$omp end parallel do
  end subroutine

  subroutine prepare(domain,constraints,resident_short)
    type(s_domain), intent(in) :: domain
    type(s_constraints), intent(in) :: constraints
    logical,optional,intent(in)::resident_short
    logical::copy_short,copy_work
    if(resident) return
    if(cells/=0) then
      if(cells/=domain%num_cell_local.or.ld/=size(domain%coord,2)) error stop 'GPU timestep shape changed'
    else
      cells=domain%num_cell_local;ld=size(domain%coord,2);workers=nthread
      allocate(host_state(3,ld,cells,7),host_cm(8,cells));host_state=0
      allocate(dx(3,ld,cells),dv(3,ld,cells),dxref(3,ld,cells),dvref(3,ld,cells),dwork(3,ld,cells))
      allocate(dfs(3,ld,cells),dfl(3,ld,cells),dmass(ld,cells),dimass(ld,cells),dcm(8,cells))
      allocate(dnatom(cells),dnwater(cells))
      allocate(dnsolute(cells),dvir(3,3,cells),host_vir(3,3,cells))
      allocate(dwater_vir(3,3,size(domain%water_list,2),cells))
      allocate(dsolute_vir(3,size(constraints%HGr_bond_list,2),constraints%connect,cells))
      allocate(dwater(size(domain%water_list,1),size(domain%water_list,2),size(domain%water_list,3)))
      allocate(dhlocal(size(constraints%HGr_local,1),size(constraints%HGr_local,2)))
      allocate(dhlist(size(constraints%HGr_bond_list,1),size(constraints%HGr_bond_list,2), &
                      size(constraints%HGr_bond_list,3),size(constraints%HGr_bond_list,4)))
      allocate(dhdist(size(constraints%HGr_bond_dist,1),size(constraints%HGr_bond_dist,2), &
                      size(constraints%HGr_bond_dist,3),size(constraints%HGr_bond_dist,4)))
    end if
    copy_short=.true.
    if(present(resident_short)) copy_short=.not.resident_short
    copy_work=.true.
    if(present(resident_short)) copy_work=.not.(resident_short.and.state_work_enabled())
    work_ready=copy_work
    call pack_state(domain,copy_short,copy_work)
    if(.not.copy_work.and.state_verify_enabled()) then
      call poison_work_c(dxref,dwork,ld,cells)
      write(*,'(a)') 'GENESIS_STATE_VERIFY: VV2 scratch inputs poisoned'
    end if
    dv=host_state(:,:,:,1);dx=host_state(:,:,:,2)
    if(copy_work) dxref=host_state(:,:,:,3)
    if(copy_work) dwork=host_state(:,:,:,4)
    dfl=host_state(:,:,:,5);dvref=host_state(:,:,:,7)
    if(copy_short) dfs=host_state(:,:,:,6)
    if(.not.maps_current) then
    dmass=domain%mass(:,1:cells);dimass=domain%inv_mass(:,1:cells)
    dnatom=domain%num_atom(1:cells);dnwater=domain%num_water(1:cells);dwater=domain%water_list
    dnsolute=constraints%No_HGr(1:cells)
    dhlocal=constraints%HGr_local;dhlist=constraints%HGr_bond_list;dhdist=constraints%HGr_bond_dist
    dcells=cells;diter=constraints%shake_iteration;dconnect=constraints%connect
    dtolerance=real(constraints%shake_tolerance,8)
    dwater_rOH=constraints%water_rOH;dwater_rHH=constraints%water_rHH
    dwater_massO=constraints%water_massO;dwater_massH=constraints%water_massH
    maps_current=.true.
    end if
    call coordinate_configure(domain,dx)
    resident=.true.
  end subroutine

  subroutine step_coords(domain)
    type(s_domain),intent(inout)::domain
    integer::cell,atom
    if(.not.resident) return
    host_state(:,:,:,2)=dx
    !$omp parallel do private(atom) schedule(static)
    do cell=1,cells
      do atom=1,domain%num_atom(cell)
        domain%coord(:,atom,cell)=host_state(:,atom,cell,2)
      end do
    end do
    !$omp end parallel do
    call coordinate_host_refreshed_c()
  end subroutine

  subroutine step_run1(dynamics,istep,dt_long,dt_short,ensemble,domain,constraints,boundary,dynvars,thermostat_ready)
    type(s_dynamics),intent(in)::dynamics
    integer,intent(in)::istep
    real(wip),intent(in)::dt_long,dt_short
    type(s_ensemble),intent(in)::ensemble
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    type(s_boundary),intent(inout)::boundary
    type(s_dynvars),intent(inout)::dynvars
    real(8)::scale(3),velocity_scale(3),bmoment(3),gr,half_short,half_long
    integer::cell,atom,bad
    logical,optional,intent(in)::thermostat_ready
    logical::long_step,prepared,short_force,download_coords
    prepared=.false.
    if(present(thermostat_ready)) prepared=thermostat_ready
    call prepare(domain,constraints)
    boundary%box_size_x_ref=boundary%box_size_x
    boundary%box_size_y_ref=boundary%box_size_y
    boundary%box_size_z_ref=boundary%box_size_z
    long_step=mod(istep-1,dynamics%elec_long_period)==0
    half_short=0.5d0*dt_short;half_long=0.5d0*dt_long
    bmoment=dynvars%barostat_momentum
    if(.not.prepared) then
    call save_reference<<<(ld*cells+255)/256,256>>>(ld,cells)
    call checked(cudaGetLastError())
    end if
    if(long_step) then
      gr=(bmoment(1)+bmoment(2)+bmoment(3))/ensemble%degree
      velocity_scale=exp(-(bmoment+gr)*half_long);dscale=velocity_scale
      if(prepared) then
        call step_thermo_scale(domain,constraints,velocity_scale)
      else
      call scale_groups<<<dim3(cells,(ld+63)/64,constraints%connect+2),64>>>(1)
      call checked(cudaGetLastError())
      call scale_groups<<<dim3(cells,(ld+63)/64,constraints%connect+2),64>>>(1)
      call checked(cudaGetLastError())
      end if
    end if
    call position_kick<<<(ld*cells+255)/256,256>>>(ld,cells,half_long,half_short,long_step)
    call checked(cudaGetLastError())
    scale=exp(bmoment*half_short);dscale=scale
    call scale_groups<<<dim3(cells,(ld+63)/64,constraints%connect+2),64>>>(2)
    call checked(cudaGetLastError())
    call position_drift<<<(ld*cells+255)/256,256>>>(ld,cells,real(dt_short,8))
    call checked(cudaGetLastError())
    call scale_groups<<<dim3(cells,(ld+63)/64,constraints%connect+2),64>>>(3)
    call checked(cudaGetLastError())
    failed=0
    call water_settle<<<dim3(cells,(size(domain%water_list,2)+63)/64,1),64>>>(real(dt_short,8))
    call checked(cudaGetLastError())
    call solute_shake<<<dim3(cells,(size(constraints%HGr_bond_list,2)+63)/64,constraints%connect),64>>>(real(dt_short,8))
    call checked(cudaGetLastError())
    call constraint_virial<<<(cells+127)/128,128>>>()
    call checked(cudaGetLastError())
    host_vir=dvir;bad=failed
    if(bad/=0) error stop 'GPU SHAKE failed to converge'
    dynvars%virial_const=0
    do cell=1,cells
      dynvars%virial_const=dynvars%virial_const+host_vir(:,:,cell)
    end do
    dynvars%virial_const=dynvars%virial_const/(dt_short*dt_short)
    scale=exp(bmoment*dt_short)
    boundary%box_size_x=scale(1)*boundary%box_size_x_ref
    boundary%box_size_y=scale(2)*boundary%box_size_y_ref
    boundary%box_size_z=scale(3)*boundary%box_size_z_ref
    boundary%cell_size_x=boundary%box_size_x/real(boundary%num_cells_x,wip)
    boundary%cell_size_y=boundary%box_size_y/real(boundary%num_cells_y,wip)
    boundary%cell_size_z=boundary%box_size_z/real(boundary%num_cells_z,wip)
    !$omp parallel do private(atom) schedule(static)
    do cell=1,cells
      do atom=1,domain%num_atom(cell)
        domain%trans_vec(:,atom,cell)=domain%trans_vec(:,atom,cell)*scale
      end do
    end do
    !$omp end parallel do
    domain%system_size=[boundary%box_size_x,boundary%box_size_y,boundary%box_size_z]
    short_force=mod(istep,dynamics%elec_long_period)/=0
    if(coordinate_enabled()) call coordinate_pack_c(scale(1),scale(2),scale(3),merge(1,0,short_force))
    download_coords=.true.
    if(coordinate_enabled().and.coordinate_resident().and.short_force) download_coords=.false.
    if(coordinate_enabled().and.coordinate_resident().and.coordinate_pme()) then
      if(dynamics%nbupdate_period>0.and.dynamics%eneout_period>0) then
        if(dynamics%baro_period==dynamics%nbupdate_period.and. &
           mod(dynamics%eneout_period,dynamics%nbupdate_period)==0) download_coords=.false.
      end if
    end if
    if(step_verify().or.coordinate_verify()==1) download_coords=.true.
    if(download_coords) call step_coords(domain)
    if(coordinate_enabled()) call coordinate_check(domain)
    if(.not.download_coords.and.coordinate_verify()==2) then
      do cell=1,cells
        domain%coord(:,1:domain%num_atom(cell),cell)=ieee_value(0.0_wip,ieee_quiet_nan)
      end do
      write(*,'(a)') 'GENESIS_COORD_VERIFY: CPU coordinates poisoned'
    end if
  end subroutine

  attributes(global) subroutine save_reference(width,count)
    integer,value::width,count
    integer::index,cell,atom
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>dnatom(cell)) return
    dxref(:,atom,cell)=dx(:,atom,cell);dvref(:,atom,cell)=dv(:,atom,cell)
  end subroutine

  attributes(global) subroutine position_kick(width,count,half_long,half_short,long_step)
    integer,value::width,count
    real(8),value::half_long,half_short
    logical,value::long_step
    integer::index,cell,atom,k
    real(8)::factor
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>dnatom(cell)) return
    if(long_step) then
      factor=half_long*dimass(atom,cell)
      dv(:,atom,cell)=dv(:,atom,cell)+factor*dfl(:,atom,cell)
    end if
    factor=half_short*dimass(atom,cell)
    dv(:,atom,cell)=dv(:,atom,cell)+factor*dfs(:,atom,cell)
  end subroutine

  attributes(global) subroutine position_drift(width,count,dt)
    integer,value::width,count
    real(8),value::dt
    integer::index,cell,atom
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>dnatom(cell)) return
    dx(:,atom,cell)=dx(:,atom,cell)+dv(:,atom,cell)*dt
  end subroutine

  attributes(global) subroutine scale_groups(mode)
    integer,value::mode
    integer::cell,entry,group,n,ids(8),h,k
    real(8)::center(3),total,value(3)
    cell=blockIdx%x;group=blockIdx%z-2
    entry=(blockIdx%y-1)*blockDim%x+threadIdx%x
    if(group==-1) then
      if(entry>dnsolute(cell)) return
      if(mode==1) dv(:,entry,cell)=dv(:,entry,cell)*dscale
      if(mode==2) dx(:,entry,cell)=dxref(:,entry,cell)*dscale
      if(mode==3) dx(:,entry,cell)=dx(:,entry,cell)*dscale
      return
    else if(group==0) then
      if(entry>dnwater(cell)) return
      n=3;ids(1:3)=dwater(:,entry,cell)
    else
      if(entry>dhlocal(group,cell)) return
      n=group+1;ids(1:n)=dhlist(1:n,entry,group,cell)
    end if
    center=0;total=0
    do h=1,n
      if(mode==1) then
        value=dv(:,ids(h),cell)
      else if(mode==2) then
        value=dxref(:,ids(h),cell)
      else
        value=dx(:,ids(h),cell)
      end if
      center=center+dmass(ids(h),cell)*value;total=total+dmass(ids(h),cell)
    end do
    center=center/total
    do h=1,n
      if(mode==1) then
        dv(:,ids(h),cell)=dv(:,ids(h),cell)+(dscale-1.0d0)*center
      else if(mode==2) then
        dx(:,ids(h),cell)=dxref(:,ids(h),cell)+(dscale-1.0d0)*center
      else
        dx(:,ids(h),cell)=dx(:,ids(h),cell)+(dscale-1.0d0)*center
      end if
    end do
  end subroutine

  attributes(global) subroutine constraint_virial()
    integer::cell,entry,j,k
    real(8)::vir(3,3)
    cell=(blockIdx%x-1)*blockDim%x+threadIdx%x
    if(cell>dcells) return
    vir=0
    do entry=1,dnwater(cell)
      vir=vir-dwater_vir(:,:,entry,cell)
    end do
    do j=1,dconnect
      do entry=1,dhlocal(j,cell)
        do k=1,3
          vir(k,k)=vir(k,k)-dsolute_vir(k,entry,j,cell)
        end do
      end do
    end do
    dvir(:,:,cell)=vir
  end subroutine

  subroutine step_sync(domain,vv1_only)
    type(s_domain), intent(inout) :: domain
    logical,optional,intent(in)::vv1_only
    integer :: cell,atom
    logical::copy_short,copy_work,copy_ref,trimmed
    integer::status
    character(16)::value
    if(.not.resident) return
    call hgr_latch_c(dx,dv,dhlist,dhlocal,dhdist,ld,cells,size(dhlist,2),size(dhlist,1),size(dhlist,3))
    call water_latch_c(dx,dv,dwater,dnwater,ld,cells,size(dwater,2))
    if(thermo_pending) call thermo_output(domain)
    copy_short=force_resident_requested();copy_work=work_ready
    trimmed=.false.
    if(present(vv1_only)) then
      call get_environment_variable('GENESIS_STATE_HANDOFF',value,status=status)
      trimmed=vv1_only.and.status==0.and.trim(value)=='1'
    end if
    copy_ref=.not.trimmed
    if(trimmed) copy_work=.false.
    host_state(:,:,:,1)=dv;host_state(:,:,:,2)=dx
    if(copy_work.and.copy_ref) host_state(:,:,:,3)=dxref
    if(copy_work) host_state(:,:,:,4)=dwork
    host_state(:,:,:,5)=dfl
    if(copy_ref) host_state(:,:,:,7)=dvref
    if(copy_short) host_state(:,:,:,6)=dfs
    !$omp parallel do private(atom) schedule(static)
    do cell=1,cells
      do atom=1,domain%num_atom(cell)
        domain%velocity(:,atom,cell)=host_state(:,atom,cell,1)
        domain%coord(:,atom,cell)=host_state(:,atom,cell,2)
        if(copy_work.and.copy_ref) domain%coord_ref(:,atom,cell)=host_state(:,atom,cell,3)
        if(copy_work) domain%velocity_full(:,atom,cell)=host_state(:,atom,cell,4)
        domain%force_long(:,atom,cell)=host_state(:,atom,cell,5)
        if(copy_short) domain%force_short(:,atom,cell)=host_state(:,atom,cell,6)
        if(copy_ref) domain%velocity_ref(:,atom,cell)=host_state(:,atom,cell,7)
      end do
    end do
    !$omp end parallel do
    call coordinate_host_refreshed_c()
    if(trimmed.and.state_verify_enabled()) then
      do cell=1,cells
        domain%coord_ref(:,1:domain%num_atom(cell),cell)=ieee_value(0.0_wip,ieee_quiet_nan)
        domain%velocity_ref(:,1:domain%num_atom(cell),cell)=ieee_value(0.0_wip,ieee_quiet_nan)
      end do
      write(*,'(a)') 'GENESIS_STATE_VERIFY: CPU VV1 references poisoned'
    end if
  end subroutine

  subroutine step_run2(dynamics,istep,dt_long,dt_short,ensemble,domain,constraints,boundary,dynvars)
    type(s_dynamics), intent(in) :: dynamics
    integer,intent(in)::istep
    real(wip),intent(in)::dt_long,dt_short
    type(s_ensemble),intent(in)::ensemble
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    type(s_boundary),intent(inout)::boundary
    type(s_dynvars),intent(in)::dynvars
    real(8)::cm(8),half_long,half_short,bmoment(3)
    logical::long_step
    integer::cell,atom,thread,status,bad
    call step_force2(dynamics,istep,domain,constraints)
    long_step=mod(istep,dynamics%elec_long_period)==0
    cm=0
    if(long_step) then
      call center_momentum<<<cells,128>>>(cells)
      call checked(cudaGetLastError())
      host_cm=dcm
      do thread=1,cells
        cm=cm+host_cm(:,thread)
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place,cm,8,mpi_real8,mpi_sum,mpi_comm_city,status)
      if(status/=MPI_SUCCESS) error stop 'GPU timestep COM reduction failed'
#endif
    end if
    half_long=0.5d0*dt_long;half_short=0.5d0*dt_short
    bmoment=dynvars%barostat_momentum
    dcenter=cm;dbmoment=bmoment
    failed=0
    call kick_derivative<<<dim3((ld*cells+255)/256,1,1),dim3(256,1,1)>>>( &
                           ld,cells,half_long,half_short,long_step)
    call checked(cudaGetLastError())
    call water_rattle<<<dim3(cells,(size(domain%water_list,2)+63)/64,1),dim3(64,1,1)>>>()
    call checked(cudaGetLastError())
    call solute_rattle<<<dim3(cells,(size(constraints%HGr_bond_list,2)+63)/64,constraints%connect),dim3(64,1,1)>>>()
    call checked(cudaGetLastError())
    call finish_derivative<<<(ld*cells+255)/256,256>>>(ld,cells)
    work_ready=.true.
    call checked(cudaGetLastError())
    bad=failed
    if(bad/=0) error stop 'GPU RATTLE failed to converge'
    boundary%box_size_x_ref=boundary%box_size_x
    boundary%box_size_y_ref=boundary%box_size_y
    boundary%box_size_z_ref=boundary%box_size_z
  end subroutine

  subroutine step_force2(dynamics,istep,domain,constraints)
    type(s_dynamics),intent(in)::dynamics
    integer,intent(in)::istep
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    integer::cell,atom,taken
    logical::long_step
    thermo_pending=.false.
    call prepare(domain,constraints,force_resident_requested())
    taken=assembly_copy(dfs,ld,cells)
    if(force_resident_requested().and.taken==0) error stop 'Resident short force is unavailable'
    long_step=mod(istep,dynamics%elec_long_period)==0
    if(taken==0.or.long_step) then
    !$omp parallel do private(atom) schedule(static)
    do cell=1,cells
      do atom=1,domain%num_atom(cell)
        if(taken==0) host_state(:,atom,cell,6)=domain%force_short(:,atom,cell)
        if(long_step) host_state(:,atom,cell,5)=domain%force_long(:,atom,cell)
      end do
    end do
    !$omp end parallel do
    end if
    if(taken==0) dfs=host_state(:,:,:,6)
    if(long_step) dfl=host_state(:,:,:,5)
  end subroutine

  attributes(global) subroutine center_momentum(threads)
    integer,value::threads
    integer::thread,cell,atom,k,stride
    real(8)::s1,s2,s3,s4,s5,s6,s7,s8,m
    real(8),shared::partial(8,128)
    thread=threadIdx%x;cell=blockIdx%x
    s1=0;s2=0;s3=0;s4=0;s5=0;s6=0;s7=0;s8=0
    do atom=thread,dnatom(cell),128
      m=dmass(atom,cell)
      s1=s1+m*dv(1,atom,cell);s2=s2+m*dv(2,atom,cell);s3=s3+m*dv(3,atom,cell)
      s4=s4+m;s5=s5+dfl(1,atom,cell);s6=s6+dfl(2,atom,cell);s7=s7+dfl(3,atom,cell);s8=s8+1
    end do
    partial(1,thread)=s1;partial(2,thread)=s2;partial(3,thread)=s3;partial(4,thread)=s4
    partial(5,thread)=s5;partial(6,thread)=s6;partial(7,thread)=s7;partial(8,thread)=s8
    call syncthreads()
    stride=64
    do while(stride>0)
      if(thread<=stride) then
        do k=1,8
          partial(k,thread)=partial(k,thread)+partial(k,thread+stride)
        end do
      end if
      call syncthreads()
      stride=stride/2
    end do
    if(thread==1) dcm(:,cell)=partial(:,1)
  end subroutine

  attributes(global) subroutine kick_derivative(width,count,half_long,half_short,long_step)
    integer,value::width,count
    real(8),value::half_long,half_short
    logical,value::long_step
    integer::index,cell,atom,k
    real(8)::factor
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>dnatom(cell)) return
    if(long_step) then
      do k=1,3
        dv(k,atom,cell)=dv(k,atom,cell)-dcenter(k)/dcenter(4)
        dfl(k,atom,cell)=dfl(k,atom,cell)-dcenter(k+4)/dcenter(8)
      end do
      factor=half_long*dimass(atom,cell)
      do k=1,3
        dv(k,atom,cell)=dv(k,atom,cell)+factor*dfl(k,atom,cell)
      end do
    end if
    factor=half_short*dimass(atom,cell)
    do k=1,3
      dv(k,atom,cell)=dv(k,atom,cell)+factor*dfs(k,atom,cell)
      dwork(k,atom,cell)=dv(k,atom,cell)+dbmoment(k)*dx(k,atom,cell)
      dxref(k,atom,cell)=dwork(k,atom,cell)
    end do
  end subroutine

  attributes(global) subroutine finish_derivative(width,count)
    integer,value::width,count
    integer::index,cell,atom,k
    real(8)::change
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>dnatom(cell)) return
    do k=1,3
      change=dwork(k,atom,cell)-dxref(k,atom,cell)
      dv(k,atom,cell)=dv(k,atom,cell)+change
    end do
  end subroutine

  subroutine state_copy(domain,array,restore)
    type(s_domain),intent(inout)::domain
    real(8),intent(inout)::array(:,:,:,:)
    logical,intent(in)::restore
    if(restore) then
      domain%velocity(:,1:ld,1:cells)=array(:,:,:,1)
      domain%coord(:,1:ld,1:cells)=array(:,:,:,2)
      domain%coord_ref(:,1:ld,1:cells)=array(:,:,:,3)
      domain%velocity_full(:,1:ld,1:cells)=array(:,:,:,4)
      domain%force_long(:,1:ld,1:cells)=array(:,:,:,5)
      domain%velocity_ref(:,1:ld,1:cells)=array(:,:,:,6)
    else
      array(:,:,:,1)=domain%velocity(:,1:ld,1:cells)
      array(:,:,:,2)=domain%coord(:,1:ld,1:cells)
      array(:,:,:,3)=domain%coord_ref(:,1:ld,1:cells)
      array(:,:,:,4)=domain%velocity_full(:,1:ld,1:cells)
      array(:,:,:,5)=domain%force_long(:,1:ld,1:cells)
      array(:,:,:,6)=domain%velocity_ref(:,1:ld,1:cells)
    end if
  end subroutine

  subroutine step_check_begin(domain,constraints,boundary)
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    type(s_boundary),intent(in)::boundary
    call prepare(domain,constraints)
    call step_sync(domain)
    if(.not.allocated(saved)) allocate(saved(3,ld,cells,6),expected(3,ld,cells,6))
    call state_copy(domain,saved,.false.)
    saved_box=[boundary%box_size_x_ref,boundary%box_size_y_ref,boundary%box_size_z_ref]
    saved_boundary=boundary;saved_system=domain%system_size
    if(.not.allocated(saved_trans)) allocate(saved_trans,mold=domain%trans_vec)
    saved_trans=domain%trans_vec
  end subroutine

  subroutine step_check_reference(domain,boundary)
    type(s_domain),intent(inout)::domain
    type(s_boundary),intent(inout)::boundary
    call state_copy(domain,expected,.false.)
    expected_box=[boundary%box_size_x_ref,boundary%box_size_y_ref,boundary%box_size_z_ref]
    expected_boundary=boundary
    call state_copy(domain,saved,.true.)
    boundary=saved_boundary;domain%system_size=saved_system;domain%trans_vec=saved_trans
    boundary%box_size_x_ref=saved_box(1);boundary%box_size_y_ref=saved_box(2);boundary%box_size_z_ref=saved_box(3)
  end subroutine

  subroutine step_check_finish(domain,constraints,boundary,position)
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    type(s_boundary),intent(in)::boundary
    logical,optional,intent(in)::position
    real(8)::ratio,error,value,reference,atol,rtol,maxratio,residual,r(3),v(3),tolerance
    integer::cell,atom,field,k,group,entry,h,first,second,worst(4)
    real(8)::worst_actual,worst_reference
    call step_sync(domain)
    call state_copy(domain,saved,.false.)
    maxratio=0
    do field=1,6
      atol=1.0d-10;rtol=1.0d-12
      if(field==5) then
        atol=1.0d-8;rtol=1.0d-10
      end if
      do cell=1,cells
        do atom=1,domain%num_atom(cell)
          do k=1,3
            value=saved(k,atom,cell,field);reference=expected(k,atom,cell,field)
            if(.not.ieee_is_finite(value).or..not.ieee_is_finite(reference)) error stop 'Non-finite step comparison'
            ratio=abs(value-reference)/(atol+rtol*abs(reference))
            if(ratio>maxratio) then
              maxratio=ratio;worst=[field,k,atom,cell];worst_actual=value;worst_reference=reference
            end if
          end do
        end do
      end do
    end do
    if(maxratio>1) then
      print *, 'GENESIS_STEP_VERIFY failed normalized state error ',maxratio
      print *, 'Worst field,component,atom,cell ',worst
      print *, 'GPU/CPU state values ',worst_actual,worst_reference
      cell=worst(4);atom=worst(3)
      do group=1,constraints%connect
        do entry=1,constraints%HGr_local(group,cell)
          if(.not.any(constraints%HGr_bond_list(1:group+1,entry,group,cell)==atom)) cycle
          first=constraints%HGr_bond_list(1,entry,group,cell)
          do h=1,group
            second=constraints%HGr_bond_list(h+1,entry,group,cell)
            r=saved(:,first,cell,2)-saved(:,second,cell,2)
            v=expected(:,first,cell,2)-expected(:,second,cell,2)
            value=constraints%HGr_bond_dist(h+1,entry,group,cell)
            print *, 'Constraint squared-distance errors GPU/CPU ',sum(r*r)-value*value,sum(v*v)-value*value
          end do
        end do
      end do
      if(.not.step_diagnose()) error stop 'CPU/GPU integration state mismatch'
    end if
    if(any(abs([boundary%box_size_x_ref,boundary%box_size_y_ref,boundary%box_size_z_ref]-expected_box)>1.0d-12)) &
      error stop 'CPU/GPU box mismatch'
    if(present(position)) then
      if(step_diagnose()) then
        write(*,'(a,es20.12)') 'GENESIS_STEP_DIAGNOSTIC: position state ratio ',maxratio
      else
        write(*,'(a,es20.12)') 'GENESIS_STEP_VERIFY: position state pass ',maxratio
      end if
      return
    end if
    residual=0
    do cell=1,cells
      do group=1,constraints%connect
        do entry=1,constraints%HGr_local(group,cell)
          first=constraints%HGr_bond_list(1,entry,group,cell)
          do h=1,group
            second=constraints%HGr_bond_list(h+1,entry,group,cell)
            r=domain%coord(:,first,cell)-domain%coord(:,second,cell)
            v=domain%velocity_full(:,first,cell)-domain%velocity_full(:,second,cell)
            residual=max(residual,abs(sum(r*v)))
          end do
        end do
      end do
      do entry=1,domain%num_water(cell)
        do first=1,2
          do second=first+1,3
            atom=domain%water_list(first,entry,cell);h=domain%water_list(second,entry,cell)
            r=domain%coord(:,atom,cell)-domain%coord(:,h,cell)
            v=domain%velocity_full(:,atom,cell)-domain%velocity_full(:,h,cell)
            residual=max(residual,abs(sum(r*v)))
          end do
        end do
      end do
    end do
    tolerance=real(constraints%shake_tolerance,8)
    if(residual>1.01d0*tolerance+1.0d-12) error stop 'GPU velocity constraint residual too large'
    if(step_diagnose()) then
      write(*,'(a,2(es20.12,1x))') 'GENESIS_STEP_DIAGNOSTIC: velocity state ratio and residual ',maxratio,residual
    else
      write(*,'(a,2(es20.12,1x))') 'GENESIS_STEP_VERIFY: state and velocity constraints pass ',maxratio,residual
    end if
  end subroutine

  subroutine step_check1_finish(domain,constraints,boundary,virial,reference)
    type(s_domain),intent(inout)::domain
    type(s_constraints),intent(in)::constraints
    type(s_boundary),intent(in)::boundary
    real(8),intent(in)::virial(3,3),reference(3,3)
    real(8)::ratio,distance,target,residual,r(3)
    integer::cell,group,entry,h,first,second,atom,other
    call step_check_finish(domain,constraints,boundary,.true.)
    if(.not.all(ieee_is_finite(virial)).or..not.all(ieee_is_finite(reference))) &
      error stop 'Non-finite position constraint virial'
    ratio=maxval(abs(virial-reference)/(1.0d-7+1.0d-10*abs(reference)))
    if(ratio>1) then
      print *, 'GENESIS_STEP_VERIFY virial error ',ratio
      print *, 'CPU constraint virial ',reference
      print *, 'GPU constraint virial ',virial
      print *, 'Constraint virial difference ',virial-reference
      if(.not.step_diagnose()) error stop 'GPU position constraint virial mismatch'
    end if
    if(any(abs([boundary%box_size_x,boundary%box_size_y,boundary%box_size_z]- &
               [expected_boundary%box_size_x,expected_boundary%box_size_y,expected_boundary%box_size_z])>1.0d-12)) &
      error stop 'GPU position box mismatch'
    residual=0
    do cell=1,cells
      do group=1,constraints%connect
        do entry=1,constraints%HGr_local(group,cell)
          first=constraints%HGr_bond_list(1,entry,group,cell)
          do h=1,group
            second=constraints%HGr_bond_list(h+1,entry,group,cell)
            target=constraints%HGr_bond_dist(h+1,entry,group,cell)
            r=domain%coord(:,first,cell)-domain%coord(:,second,cell)
            residual=max(residual,abs(sum(r*r)-target*target)/(2*constraints%shake_tolerance*target+1.0d-12))
          end do
        end do
      end do
      do entry=1,domain%num_water(cell)
        do first=1,2
          do second=first+1,3
            atom=domain%water_list(first,entry,cell);other=domain%water_list(second,entry,cell)
            target=constraints%water_rOH
            if(first==2) target=constraints%water_rHH
            r=domain%coord(:,atom,cell)-domain%coord(:,other,cell)
            residual=max(residual,abs(sum(r*r)-target*target)/(2*constraints%shake_tolerance*target+1.0d-12))
          end do
        end do
      end do
    end do
    if(residual>1.01d0) error stop 'GPU position constraint residual too large'
    if(step_diagnose()) then
      write(*,'(a,2(es20.12,1x))') 'GENESIS_STEP_DIAGNOSTIC: position residual and virial ratio ',residual,ratio
    else
      write(*,'(a,2(es20.12,1x))') 'GENESIS_STEP_VERIFY: position constraints and virial pass ',residual,ratio
    end if
  end subroutine

  attributes(global) subroutine water_rattle()
    ! local variables
    real(dp)                  :: rHH, rOH
    real(dp)                  :: massO, massH, massOH, massHH, mass2OH
    real(dp)                  :: rab(1:3), rbc(1:3), rca(1:3)
    real(dp)                  :: vab(1:3), vbc(1:3), vca(1:3)
    real(dp)                  :: unit_ab(1:3), unit_bc(1:3), unit_ca(1:3)
    real(dp)                  :: length_ab, length_bc, length_ca
    real(dp)                  :: vab0, vbc0, vca0, cosA, cosB, cosC
    real(dp)                  :: MaBC, MbCA, McAB, T_AB, T_BC, T_CA, Det
    integer                   :: i, ix, iatom(1:3), k
    integer                   :: id, omp_get_thread_num


    rOH        = dwater_rOH
    rHH        = dwater_rHH
    massO      = dwater_massO
    massH      = dwater_massH

    massOH     = massO + massH
    massHH     = massH + massH
    mass2OH    = 2.0_dp * massOH


    i=blockIdx%x
    ix=(blockIdx%y-1)*blockDim%x+threadIdx%x
    if(i>dcells) return
    if(ix>dnwater(i)) return

      iatom(1:3) = dwater(1:3,ix,i)

      rab(1:3) = dx(1:3,iatom(2),i) - dx(1:3,iatom(1),i)
      rbc(1:3) = dx(1:3,iatom(3),i) - dx(1:3,iatom(2),i)
      rca(1:3) = dx(1:3,iatom(1),i)  - dx(1:3,iatom(3),i)

      vab(1:3) = dwork(1:3,iatom(2),i) - dwork(1:3,iatom(1),i)
      vbc(1:3) = dwork(1:3,iatom(3),i) - dwork(1:3,iatom(2),i)
      vca(1:3) = dwork(1:3,iatom(1),i)  - dwork(1:3,iatom(3),i)

      length_ab = sqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)
      length_bc = sqrt(rbc(1)**2 + rbc(2)**2 + rbc(3)**2)
      length_ca = sqrt(rca(1)**2 + rca(2)**2 + rca(3)**2)

      unit_ab(1:3) = rab(1:3)/length_ab
      unit_bc(1:3) = rbc(1:3)/length_bc
      unit_ca(1:3) = rca(1:3)/length_ca

      vab0 = 0.0_dp
      vbc0 = 0.0_dp
      vca0 = 0.0_dp

      do k = 1, 3
        vab0 = vab0 + vab(k)*unit_ab(k)
        vbc0 = vbc0 + vbc(k)*unit_bc(k)
        vca0 = vca0 + vca(k)*unit_ca(k)
      end do

      cosA = 0.0_dp
      cosB = 0.0_dp
      cosC = 0.0_dp

      do k = 1, 3
        cosA = cosA - unit_ab(k)*unit_ca(k)
        cosB = cosB - unit_bc(k)*unit_ab(k)
        cosC = cosC - unit_ca(k)*unit_bc(k)
      end do

      MbCA = massH*cosC*cosA - massOH*cosB
      MaBC = massO*cosB*cosC - massHH*cosA
      McAB = massH*cosA*cosB - massOH*cosC

      T_AB = vab0*(mass2OH-massO*cosC**2) + vbc0*MbCA + vca0*MaBC
      T_AB = T_AB * massO
      T_BC = vbc0*(massOH**2-(massH*cosA)**2)
      T_BC = T_BC + vca0*massO*McAB + vab0*massO*MbCA
      T_CA = vca0*(mass2OH-massO*cosB**2) + vab0*MaBC + vbc0*McAB
      T_CA = T_CA * massO
      Det  = 2.0_dp*massOH**2 + 2.0_dp*massO*massH*cosA*cosB*cosC
      Det  = Det - 2.0_dp*massH**2*cosA**2 - massO*massOH*(cosB**2+cosC**2)
      Det  = Det / massH

      T_AB = T_AB / Det
      T_BC = T_BC / Det
      T_CA = T_CA / Det

      dwork(1:3,iatom(1),i)  = dwork(1:3,iatom(1),i) +  &
                       real((T_AB*unit_ab(1:3)-T_CA*unit_ca(1:3))/massO,wip)
      dwork(1:3,iatom(2),i) = dwork(1:3,iatom(2),i) + &
                       real((T_BC*unit_bc(1:3)-T_AB*unit_ab(1:3))/massH,wip)
      dwork(1:3,iatom(3),i) = dwork(1:3,iatom(3),i) + &
                       real((T_CA*unit_ca(1:3)-T_BC*unit_bc(1:3))/massH,wip)

  end subroutine
  attributes(global) subroutine solute_rattle()
    ! local variables
    real(dp)                     :: tolerance, imass1, imass2
    real(dp)                     :: x12, y12, z12
    real(dp)                     :: vel_x12, vel_y12, vel_z12
    real(dp)                     :: r, dot_d12_v12
    real(dp)                     :: factor, g12, g12m1, g12m2
    real(dp)                     :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: iteration, omp_get_thread_num
    integer                      :: iatm(1:8)
    logical                      :: rattle_end


    icel=blockIdx%x;j=blockIdx%z
    k=(blockIdx%y-1)*blockDim%x+threadIdx%x
    if(icel>dcells.or.j>dconnect) return
    if(k>dhlocal(j,icel)) return
    iteration=diter;tolerance=dtolerance

          iatm(1:j+1) = dhlist(1:j+1,k,j,icel)
          imass1 = 1.0_dp / dmass(iatm(1),icel)

          coord_dtmp(1:3,1:j+1) = &
            real(dx(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(dwork  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            rattle_end = .true.

            do ih = 1, j

              imass2 = 1.0_dp / dble(dmass(iatm(ih+1),icel))

              r = real(dhdist(ih+1,k,j,icel),dp)

              x12     = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12     = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12     = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              vel_x12 = vel_dtmp(1,1)   - vel_dtmp(1,ih+1)
              vel_y12 = vel_dtmp(2,1)   - vel_dtmp(2,ih+1)
              vel_z12 = vel_dtmp(3,1)   - vel_dtmp(3,ih+1)

              dot_d12_v12 = x12*vel_x12 + y12*vel_y12 + z12*vel_z12

              if (abs(dot_d12_v12) >= tolerance) then

                rattle_end = .false.

                g12 = dot_d12_v12 / ( (imass1+imass2)*r**2)
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2

                vel_dtmp(1,1)    = vel_dtmp(1,1)    - g12m1 * x12
                vel_dtmp(2,1)    = vel_dtmp(2,1)    - g12m1 * y12
                vel_dtmp(3,1)    = vel_dtmp(3,1)    - g12m1 * z12
                vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + g12m2 * x12
                vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + g12m2 * y12
                vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + g12m2 * z12

              end if
            end do

            if (rattle_end) exit

          end do

          dx(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),wip)

          dwork  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),wip)

          if (.not.rattle_end) id=atomicadd(failed,1)

  end subroutine
  attributes(global) subroutine water_settle(dt)
    real(8),value::dt
    ! local variables
    real(dp)                 :: rHH, rOH, massO, massH
    real(dp)                 :: mass(1:3), mass_H2O, ra, rb, rc, inv_ra, inv_dt
    real(dp)                 :: x0(1:3,1:3), x1(1:3,1:3), x3(1:3,1:3)
    real(dp)                 :: delt(1:3,1:3), rf(1:3,1:3,1:3)
    real(dp)                 :: com0(1:3), com1(1:3)
    real(dp)                 :: oh1(1:3), oh2(1:3)
    real(dp)                 :: Xaxis(1:3), Yaxis(1:3), Zaxis(1:3)
    real(dp)                 :: Xaxis2(1:3), Yaxis2(1:3), Zaxis2(1:3)
    real(dp)                 :: mtrx(1:3,1:3)
    real(dp)                 :: xp0(1:3,1:3), xp1(1:3,1:3), xp2(1:3,1:3)
    real(dp)                 :: xp3(1:3,1:3)
    real(dp)                 :: dxp21(1:2), dxp23(1:2), dxp31(1:2)
    real(dp)                 :: sin_phi, cos_phi, sin_psi, cos_psi, tmp
    real(dp)                 :: rb_cosphi, rb_sinphi
    real(dp)                 :: rc_sinpsi_sinphi, rc_sinpsi_cosphi
    real(dp)                 :: alpha, beta, gamma, al2bt2, sin_theta, cos_theta
    real(dp)                 :: viri_local(1:3,1:3)
    integer                  :: i, ix, iatom(1:3), i1, j1, k1
    integer                  :: id, omp_get_thread_num

    rOH        =  dwater_rOH
    rHH        =  dwater_rHH
    massO      =  dwater_massO
    massH      =  dwater_massH

    mass(1)    =  massO
    mass(2)    =  massH
    mass(3)    =  massH
    mass_H2O   =  massO + 2.0_dp*massH
    rc         =  0.5_dp * rHH 
    ra         =  2.0_dp * sqrt(rOH*rOH-rc*rc) * massH / mass_H2O
    rb         =  sqrt(rOH*rOH-rc*rc) - ra
    inv_ra     =  1.0_dp / ra
    inv_dt     =  1.0_dp / dt


    i=blockIdx%x;ix=(blockIdx%y-1)*blockDim%x+threadIdx%x
    if(ix>dnwater(i)) return
    viri_local=0

        iatom(1:3) = dwater(1:3,ix,i)
        com0(1:3)   = 0.0_dp
        com1(1:3)  = 0.0_dp
        do k1 = 1, 3
          x0(1:3,k1) = dxref(1:3,iatom(k1),i)
          x1(1:3,k1) = dx    (1:3,iatom(k1),i)
          com0(1:3)  = com0(1:3) + x0(1:3,k1)*mass(k1)
          com1(1:3)  = com1(1:3) + x1(1:3,k1)*mass(k1)
        end do
        com0(1:3) = com0(1:3) / mass_H2O
        com1(1:3) = com1(1:3) / mass_H2O
  
        do k1 = 1, 3
          x0(1:3,k1) = x0(1:3,k1) - com0(1:3)
          x1(1:3,k1) = x1(1:3,k1) - com1(1:3)
        end do
  
        oh1(1:3) = x0(1:3,2) - x0(1:3,1)
        oh2(1:3) = x0(1:3,3) - x0(1:3,1)
  
        Zaxis(1) = oh1(2)*oh2(3) - oh1(3)*oh2(2)
        Zaxis(2) = oh1(3)*oh2(1) - oh1(1)*oh2(3)
        Zaxis(3) = oh1(1)*oh2(2) - oh1(2)*oh2(1)
        Xaxis(1) = x1(2,1)*Zaxis(3) - x1(3,1)*Zaxis(2)
        Xaxis(2) = x1(3,1)*Zaxis(1) - x1(1,1)*Zaxis(3)
        Xaxis(3) = x1(1,1)*Zaxis(2) - x1(2,1)*Zaxis(1)
        Yaxis(1) = Zaxis(2)*Xaxis(3) - Zaxis(3)*Xaxis(2)
        Yaxis(2) = Zaxis(3)*Xaxis(1) - Zaxis(1)*Xaxis(3)
        Yaxis(3) = Zaxis(1)*Xaxis(2) - Zaxis(2)*Xaxis(1)
  
        Xaxis2(1:3) = Xaxis(1:3) * Xaxis(1:3)
        Yaxis2(1:3) = Yaxis(1:3) * Yaxis(1:3)
        Zaxis2(1:3) = Zaxis(1:3) * Zaxis(1:3)
        mtrx(1:3,1) = Xaxis(1:3) / sqrt(Xaxis2(1)+Xaxis2(2)+Xaxis2(3))
        mtrx(1:3,2) = Yaxis(1:3) / sqrt(Yaxis2(1)+Yaxis2(2)+Yaxis2(3))
        mtrx(1:3,3) = Zaxis(1:3) / sqrt(Zaxis2(1)+Zaxis2(2)+Zaxis2(3))

        xp0(1:3,1:3) = 0.0_dp
        xp1(1:3,1:3) = 0.0_dp
        do i1 = 1, 3
          do k1 = 1, 3
            do j1 = 1, 3
              xp0(k1,i1) = xp0(k1,i1) + mtrx(j1,i1)*x0(j1,k1)
              xp1(k1,i1) = xp1(k1,i1) + mtrx(j1,i1)*x1(j1,k1)
            end do
          end do
        end do

        sin_phi = xp1(1,3)*inv_ra
        tmp = 1.0_dp - sin_phi*sin_phi
        if (tmp > 0.0_dp) then
          cos_phi = sqrt(tmp)
        else
          cos_phi = 0.0_dp
        end if

        sin_psi = (xp1(2,3) - xp1(3,3))/(rHH*cos_phi)
        tmp = 1.0_dp - sin_psi*sin_psi
        if (tmp > 0.0_dp) then
          cos_psi = sqrt(tmp)
        else
          cos_psi = 0.0_dp
        end if

        rb_cosphi = rb * cos_phi
        rb_sinphi = rb * sin_phi
        rc_sinpsi_sinphi = rc * sin_psi * sin_phi
        rc_sinpsi_cosphi = rc * sin_psi * cos_phi
  
        xp2(1,2) =   ra * cos_phi
        xp2(1,3) =   ra * sin_phi
        xp2(2,1) = - rc * cos_psi
        xp2(2,2) = - rb_cosphi - rc_sinpsi_sinphi
        xp2(2,3) = - rb_sinphi + rc_sinpsi_cosphi
        xp2(3,1) = - xp2(2,1)
        xp2(3,2) = - rb_cosphi + rc_sinpsi_sinphi
        xp2(3,3) = - rb_sinphi - rc_sinpsi_cosphi
  
        dxp21(1:2) = xp0(2,1:2) - xp0(1,1:2)
        dxp23(1:2) = xp0(2,1:2) - xp0(3,1:2)
        dxp31(1:2) = xp0(3,1:2) - xp0(1,1:2)
        alpha =   xp2(2,1)*dxp23(1) + dxp21(2)*xp2(2,2) + dxp31(2)*xp2(3,2)
        beta  = - xp2(2,1)*dxp23(2) + dxp21(1)*xp2(2,2) + dxp31(1)*xp2(3,2)
        gamma =   dxp21(1) * xp1(2,2) - xp1(2,1) * dxp21(2) &
                + dxp31(1) * xp1(3,2) - xp1(3,1) * dxp31(2)
  
        al2bt2 = alpha*alpha + beta*beta
        sin_theta = (alpha*gamma - beta*sqrt(al2bt2 - gamma*gamma))/al2bt2
        tmp = 1.0_dp - sin_theta*sin_theta
        if (tmp > 0.0_dp) then
          cos_theta = sqrt(tmp)
        else
          cos_theta = 0.0_dp
        end if

        xp3(1,1) = - xp2(1,2)*sin_theta
        xp3(1,2) =   xp2(1,2)*cos_theta
        xp3(1,3) =   xp2(1,3)
        xp3(2,1) =   xp2(2,1)*cos_theta - xp2(2,2)*sin_theta
        xp3(2,2) =   xp2(2,1)*sin_theta + xp2(2,2)*cos_theta
        xp3(2,3) =   xp2(2,3)
        xp3(3,1) =   xp2(3,1)*cos_theta - xp2(3,2)*sin_theta
        xp3(3,2) =   xp2(3,1)*sin_theta + xp2(3,2)*cos_theta
        xp3(3,3) =   xp2(3,3)
  
        x3(1:3,1:3) = 0.0_dp
        do k1 = 1, 3
          do i1 = 1, 3
            do j1 = 1, 3
              x3(i1,k1) = x3(i1,k1) + mtrx(i1,j1)*xp3(k1,j1)
            end do
          end do
          dx(1:3,iatom(k1),i) = real(x3(1:3,k1) + com1(1:3),wip)
        end do
  
        do k1 = 1, 3
          delt(1:3,k1) = x3(1:3,k1) - x1(1:3,k1)
          do i1 = 1, 3
            do j1 = 1, 3
              rf(i1,j1,k1) = - dxref(i1,iatom(k1),i)*delt(j1,k1)*mass(k1)
            end do
          end do
        end do
  
        do i1 = 1, 3
          do j1 = 1, 3
            viri_local(1:3,j1) = viri_local(1:3,j1) + rf(1:3,j1,i1)
          end do
        end do
  
        if (.true.) then
          do k1 = 1, 3
            dv(1:3,iatom(k1),i) = dv(1:3,iatom(k1),i) &
                                 + real(delt(1:3,k1)*inv_dt,wip)
          end do
        end if


    dwater_vir(:,:,ix,i)=viri_local
  end subroutine
  attributes(global) subroutine solute_shake(dt)
    real(8),value::dt
    real(8)::bond_vector(3,7),shake_force(7)
    ! local variables
    real(dp )                    :: tolerance, imass1, imass2
    real(dp )                    :: x12, y12, z12
    real(dp )                    :: x12_old, y12_old, z12_old
    real(dp )                    :: dist2, r
    real(dp )                    :: factor, g12, g12m1, g12m2, v12m1, v12m2
    real(dp )                    :: viri(1:3), fx, fy, fz
    real(dp )                    :: coord_dtmp(1:3,1:8), vel_dtmp(1:3,1:8)
    logical                      :: shake_end
    integer                      :: icel, i, j, k, ih, connect, id
    integer                      :: atm1, atm2, iatm(1:8)
    integer                      :: iteration, omp_get_thread_num


    icel=blockIdx%x;j=blockIdx%z;k=(blockIdx%y-1)*blockDim%x+threadIdx%x
    if(k>dhlocal(j,icel)) return
    iteration=diter;tolerance=dtolerance
    atm1=dhlist(1,k,j,icel);shake_force=0
    do ih=1,j
      atm2=dhlist(ih+1,k,j,icel)
      bond_vector(:,ih)=dxref(:,atm1,icel)-dxref(:,atm2,icel)
    end do

          iatm(1:j+1) = dhlist(1:j+1,k,j,icel)
          imass1 = 1.0_dp / real(dmass(iatm(1),icel),dp)

          coord_dtmp(1:3,1:j+1) = &
            real(dx(1:3,iatm(1:j+1),icel),dp)

          vel_dtmp  (1:3,1:j+1) = &
            real(dv  (1:3,iatm(1:j+1),icel),dp)

          do i = 1, iteration

            shake_end = .true.
            do ih = 1, j

              imass2 = 1.0_dp / real(dmass(iatm(ih+1),icel),dp)

              r = real(dhdist(ih+1,k,j,icel),dp)

              x12 = coord_dtmp(1,1) - coord_dtmp(1,ih+1)
              y12 = coord_dtmp(2,1) - coord_dtmp(2,ih+1)
              z12 = coord_dtmp(3,1) - coord_dtmp(3,ih+1)
              dist2 = x12**2 + y12**2 + z12**2

              if (abs(dist2 - r*r) >= 2.0_dp*tolerance*r) then 

                shake_end = .false.
                x12_old = bond_vector(1,ih)
                y12_old = bond_vector(2,ih)
                z12_old = bond_vector(3,ih)

                factor = (x12*x12_old + y12*y12_old + z12*z12_old)* &
                         (imass1 + imass2)
                g12    = 0.5_dp*(dist2 - r*r)/factor
                g12m1  = g12 * imass1
                g12m2  = g12 * imass2
                v12m1  = g12m1 / real(dt,dp)
                v12m2  = g12m2 / real(dt,dp)

                shake_force(ih) = &
                     shake_force(ih) + real(g12,wip)

                coord_dtmp(1,1)    = coord_dtmp(1,1)    - g12m1 * x12_old
                coord_dtmp(2,1)    = coord_dtmp(2,1)    - g12m1 * y12_old
                coord_dtmp(3,1)    = coord_dtmp(3,1)    - g12m1 * z12_old
                coord_dtmp(1,ih+1) = coord_dtmp(1,ih+1) + g12m2 * x12_old
                coord_dtmp(2,ih+1) = coord_dtmp(2,ih+1) + g12m2 * y12_old
                coord_dtmp(3,ih+1) = coord_dtmp(3,ih+1) + g12m2 * z12_old

                if (.true.) then
                  vel_dtmp(1,1)    = vel_dtmp(1,1)    - v12m1 * x12_old
                  vel_dtmp(2,1)    = vel_dtmp(2,1)    - v12m1 * y12_old
                  vel_dtmp(3,1)    = vel_dtmp(3,1)    - v12m1 * z12_old
                  vel_dtmp(1,ih+1) = vel_dtmp(1,ih+1) + v12m2 * x12_old
                  vel_dtmp(2,ih+1) = vel_dtmp(2,ih+1) + v12m2 * y12_old
                  vel_dtmp(3,ih+1) = vel_dtmp(3,ih+1) + v12m2 * z12_old
                end if 

              end if
            end do

            if (shake_end) exit

          end do

          if (.not.shake_end) id=atomicadd(failed,1)

          dx(1:3,iatm(1:j+1),icel) = &
               real(coord_dtmp(1:3,1:j+1),wip)

          dv  (1:3,iatm(1:j+1),icel) = &
               real(vel_dtmp  (1:3,1:j+1),wip)


    viri=0
    do ih=1,j
      viri=viri+shake_force(ih)*bond_vector(:,ih)*bond_vector(:,ih)
    end do
    dsolute_vir(:,k,j,icel)=viri
  end subroutine

end module
