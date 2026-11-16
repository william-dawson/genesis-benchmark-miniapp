#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module translation_rebuild_mod
  use incoming_gpu_mod,only:incoming_boundary_suspend
  use iso_c_binding,only:c_int,c_ptr,c_loc,c_size_t
  use iso_fortran_env,only:int8
  use translation_residency_mod,only:translation_mode
  use constants_mod
  use mpi_parallel_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_communicate_mod,only:s_comm,update_cell_boundary
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_finite
  implicit none
  private
  public::translation_rebuild_begin,translation_rebuild_finish,translation_output_skip
  real(wp),allocatable,target,save::observed(:,:,:)
  interface
    integer(c_int) function compare_c(a,b,bytes) bind(C,name='memcmp')
      import c_int,c_ptr,c_size_t
      type(c_ptr),value::a,b
      integer(c_size_t),value::bytes
    end function
    integer(c_int) function mode_c() bind(C,name='genesis_translation_rebuild_mode')
      import c_int
    end function
    integer(c_int) function begin_c() bind(C,name='genesis_translation_rebuild_begin')
      import c_int
    end function
    integer(c_int) function pending_c() bind(C,name='genesis_translation_rebuilding')
      import c_int
    end function
    subroutine poison_c() bind(C,name='genesis_translation_poisoned')
    end subroutine
    subroutine end_c() bind(C,name='genesis_translation_rebuild_end')
    end subroutine
  end interface
contains
  logical function translation_output_skip(domain,active)
    type(s_domain),intent(inout)::domain
    logical,intent(in)::active
    character(24),parameter::diagnostics(3)=[character(24):: &
      'GENESIS_FORCE_VERIFY','GENESIS_COORD_VERIFY','GENESIS_STEP_VERIFY']
    character(16)::value
    integer::check,status
    translation_output_skip=.false.
    if(.not.active.or.translation_mode()==0.or.mode_c()==0) return
    do check=1,size(diagnostics)
      call get_environment_variable(trim(diagnostics(check)),value,status=status)
      if(status==0.and.len_trim(value)>0.and.trim(value)/='0') return
    end do
    translation_output_skip=.true.
    if(mode_c()==2) then
      call poison_c()
      domain%trans_vec=ieee_value(0.0_wp,ieee_quiet_nan)
      write(*,'(a)') 'GENESIS_TRANSLATION_REBUILD: output mirror poisoned'
    end if
  end function

  subroutine translation_rebuild_begin(domain,comm)
    type(s_domain),intent(inout)::domain
    type(s_comm),intent(in)::comm
    if(mode_c()==0) return
    if(nproc_country/=1.or.domain%num_cell_boundary/=0.or.any(comm%num_cell/=0)) return
    if(begin_c()==0) return
    if(mode_c()==2) then
      domain%trans_vec=ieee_value(0.0_wp,ieee_quiet_nan)
      write(*,'(a)') 'GENESIS_TRANSLATION_REBUILD: old mirror poisoned'
    end if
  end subroutine

  subroutine translation_rebuild_finish(domain,comm,boundary,constraints)
    type(s_domain),target,intent(inout)::domain
    type(s_comm),intent(inout)::comm
    type(s_boundary),intent(inout)::boundary
    type(s_constraints),intent(inout)::constraints
    integer::cell,atom
    if(pending_c()==0) return
    if(mode_c()==2) then
      if(.not.allocated(observed)) allocate(observed,mold=domain%trans_vec)
      if(any(shape(observed)/=shape(domain%trans_vec))) error stop 'Translation rebuild shape changed'
      observed=domain%trans_vec
      do cell=1,domain%num_cell_local
        do atom=1,domain%num_atom(cell)
          if(.not.all(ieee_is_finite(observed(:,atom,cell)))) error stop 'Translation rebuild left a poisoned active value'
        end do
      end do
      ! No halos: the boundary routine reads unchanged coordinates and overwrites translations.
      domain%trans_vec=0.0_wp
      call incoming_boundary_suspend(.true.)
      call update_cell_boundary(domain,comm,boundary,constraints)
      call incoming_boundary_suspend(.false.)
      do cell=1,domain%num_cell_local
        do atom=1,domain%num_atom(cell)
          if(compare_c(c_loc(domain%trans_vec(1,atom,cell)),c_loc(observed(1,atom,cell)), &
                       int(3*storage_size(1.0_wp)/8,c_size_t))/=0) then
            write(*,'(a,3i10,6es24.16)') 'GENESIS_TRANSLATION_REBUILD: CPU/observed ', &
              cell,atom,domain%id_l2g(atom,cell),domain%trans_vec(:,atom,cell),observed(:,atom,cell)
            write(*,'(a,*(1x,z2.2))') 'CPU bytes ',transfer(domain%trans_vec(:,atom,cell),[0_int8])
            write(*,'(a,*(1x,z2.2))') 'Observed bytes ',transfer(observed(:,atom,cell),[0_int8])
            error stop 'Translation rebuild depends on retired mirror'
          end if
        end do
      end do
      domain%trans_vec=observed
      write(*,'(a)') 'GENESIS_TRANSLATION_REBUILD: exact overwrite pass'
    end if
    call end_c()
  end subroutine
end module
