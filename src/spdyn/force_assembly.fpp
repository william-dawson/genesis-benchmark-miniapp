#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module force_assembly_mod
  use iso_c_binding
  use constants_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  public :: force_assembly_enabled,force_assembly_run
  public :: force_resident_set,force_resident_active,force_resident_requested
  logical,save::resident_step=.false.
  real(wip),allocatable :: reference(:,:,:),observed(:,:,:)
  interface
    subroutine assemble_c(cpu,fl,output,natom,start,width,cells,capacity,outcells,workers,atoms,merge,mode,download) &
        bind(C,name='genesis_force_assemble')
      import c_double,c_int
      real(c_double),intent(in)::cpu(*),fl(*)
      real(c_double),intent(out)::output(*)
      integer(c_int),intent(in)::natom(*),start(*)
      integer(c_int),value::width,cells,capacity,outcells,workers,atoms,merge,mode,download
    end subroutine
  end interface
contains
  logical function force_resident_requested()
    integer,save::selected=-1
    integer::status,code
    character(16)::value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_FORCE_RESIDENT',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0.or.selected<0.or.selected>1) error stop 'GENESIS_FORCE_RESIDENT must be 0 or 1'
      end if
    end if
    force_resident_requested=selected==1
  end function

  subroutine force_resident_set(eligible)
    logical,intent(in)::eligible
    resident_step=eligible.and.force_resident_requested().and.force_assembly_enabled()
  end subroutine

  logical function force_resident_active()
    force_resident_active=resident_step
  end function
  integer function assembly_mode() result(selected)
    integer,save::cached=-1
    character(16)::value
    integer::status,code
    if(cached<0) then
      cached=0
      call get_environment_variable('GENESIS_FORCE_ASSEMBLY',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) cached
        if(code/=0.or.cached<0.or.cached>2) error stop 'GENESIS_FORCE_ASSEMBLY must be 0, 1 or 2'
      end if
    end if
    selected=cached
  end function

  logical function force_assembly_enabled()
    force_assembly_enabled=assembly_mode()/=0
    if(force_assembly_enabled) then
      if(nproc_country/=1) error stop 'GPU force assembly requires one MPI rank'
#if defined(_SINGLE) || defined(_MIXED)
      error stop 'GPU force assembly requires double precision'
#endif
    end if
  end function

  subroutine force_assembly_run(domain,slabs,pbc,fl,merge,output)
    type(s_domain),intent(in)::domain
    real(wp),contiguous,intent(in)::slabs(:,:,:,:)
    real(wp),contiguous,intent(in)::pbc(:,:,:,:)
    real(wip),contiguous,intent(in)::fl(:,:,:)
    logical,intent(in)::merge
    real(wip),contiguous,intent(inout)::output(:,:,:)
    integer::download,use_long,cell,atom,worker,start,atoms,status
    real(8)::x,y,z,ratio,maxratio
    character(16)::value
    logical::verify
    if(domain%num_cell_boundary/=0) error stop 'GPU force assembly halo unsupported'
    if(size(slabs,2)/=size(output,2)) error stop 'GPU force assembly pitch mismatch'
    use_long=0
    if(merge) use_long=1
    call get_environment_variable('GENESIS_FORCE_VERIFY',value,status=status)
    verify=status==0.and.trim(value)=='1'
    atoms=domain%num_atom_domain
    if(verify) then
      if(.not.allocated(reference)) allocate(reference(size(output,1),size(output,2),size(output,3)))
      !$omp parallel do private(atom,worker,start,x,y,z) schedule(static)
      do cell=1,domain%num_cell_local
        start=domain%start_atom(cell)
        do atom=1,domain%num_atom(cell)
          x=slabs(1,atom,cell,1);y=slabs(2,atom,cell,1);z=slabs(3,atom,cell,1)
          if(merge) then
            x=x+fl(1,atom,cell);y=y+fl(2,atom,cell);z=z+fl(3,atom,cell)
          end if
          x=x+pbc(start+atom,1,1,1);y=y+pbc(atoms+start+atom,1,1,1);z=z+pbc(2*atoms+start+atom,1,1,1)
          do worker=2,size(slabs,4)
            x=x+slabs(1,atom,cell,worker);y=y+slabs(2,atom,cell,worker);z=z+slabs(3,atom,cell,worker)
          end do
          reference(:,atom,cell)=[x,y,z]
        end do
      end do
      !$omp end parallel do
    end if
    download=1
    if(resident_step) download=0
    if(verify) then
      if(.not.allocated(observed)) allocate(observed,mold=output)
    call assemble_c(slabs,fl,observed,domain%num_atom,domain%start_atom,size(output,2), &
        domain%num_cell_local,size(slabs,3),size(output,3),size(slabs,4),domain%num_atom_domain,use_long,assembly_mode(),1)
      if(resident_step) then
        output=ieee_value(0.0_wip,ieee_quiet_nan)
      else
        output=observed
      end if
    else
    call assemble_c(slabs,fl,output,domain%num_atom,domain%start_atom,size(output,2), &
        domain%num_cell_local,size(slabs,3),size(output,3),size(slabs,4),domain%num_atom_domain,use_long,assembly_mode(),download)
    end if
    if(verify) then
      maxratio=0
      do cell=1,domain%num_cell_local
        do atom=1,domain%num_atom(cell)
          if(.not.all(ieee_is_finite(observed(:,atom,cell)))) error stop 'GPU assembly nonfinite force'
          ratio=maxval(abs(observed(:,atom,cell)-reference(:,atom,cell))/(1.0d-8+1.0d-10*abs(reference(:,atom,cell))))
          maxratio=max(maxratio,ratio)
        end do
      end do
      if(maxratio>1) error stop 'GPU force assembly disagrees with CPU sum'
      write(*,'(a,es20.12)') 'GENESIS_FORCE_VERIFY: force assembly pass ',maxratio
    end if
  end subroutine
end module
