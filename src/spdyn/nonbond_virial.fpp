#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module nonbond_virial_mod
  use iso_c_binding,only:c_int,c_double
  use sp_domain_str_mod
  use force_assembly_mod,only:force_assembly_enabled
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::nonbond_virial_prepare,nonbond_virial_begin,nonbond_virial_end,nonbond_energy_cpu_required
  real(c_double),allocatable,save::candidate(:,:,:)
  integer,save::active_mode=0
  interface
    subroutine request_c(eligible,atoms,needed) bind(C,name="genesis_nonbond_virial_request")
      import c_int
      integer(c_int),value::eligible,atoms,needed
    end subroutine
    integer(c_int) function get_c(values) bind(C,name="genesis_nonbond_virial_get")
      import c_int,c_double
      real(c_double),intent(out)::values(3)
    end function
    integer(c_int) function required_c() bind(C,name="genesis_nonbond_virial_host_required")
      import c_int
    end function
    subroutine poison_c(coordinates,force) bind(C,name="genesis_nonbond_virial_poison")
      import c_double
      real(c_double),intent(inout)::coordinates(*),force(*)
    end subroutine
  end interface
contains
  subroutine nonbond_virial_prepare(domain,needed)
    type(s_domain),intent(in)::domain
    logical,intent(in)::needed
    logical::eligible
    eligible=force_assembly_enabled().and.domain%num_cell_boundary==0
    call request_c(merge(1,0,eligible),domain%num_atom_domain,merge(1,0,needed))
  end subroutine
  logical function nonbond_energy_cpu_required()
    nonbond_energy_cpu_required=required_c()/=0
  end function
  logical function nonbond_virial_begin(virial) result(handled)
    real(c_double),intent(inout)::virial(:,:,:)
    real(c_double)::values(3)
    integer::axis
    active_mode=get_c(values)
    handled=active_mode==1
    if(active_mode==0)return
    if(active_mode==2)then
      if(.not.allocated(candidate))allocate(candidate,mold=virial)
      candidate=virial
      do axis=1,3
        candidate(axis,axis,1)=candidate(axis,axis,1)+values(axis)
      end do
    else
      do axis=1,3
        virial(axis,axis,1)=virial(axis,axis,1)+values(axis)
      end do
    end if
  end function
  subroutine nonbond_virial_end(virial,coordinates,force)
    real(c_double),intent(inout)::virial(:,:,:),coordinates(:,:,:),force(:,:,:,:)
    real(c_double)::actual(3,3),reference(3,3),ratio
    if(active_mode/=2)return
    actual=sum(candidate,dim=3);reference=sum(virial,dim=3)
    if(.not.all(ieee_is_finite(actual)).or..not.all(ieee_is_finite(reference))) &
      error stop 'Nonfinite nonbond virial comparison'
    ratio=maxval(abs(actual-reference)/(1.0d-7+1.0d-10*abs(reference)))
    if(ratio>1)error stop 'GPU nonbond virial disagrees with CPU sum'
    virial=candidate
    call poison_c(coordinates(1,1,1),force(1,1,1,1))
    write(*,'(a,es20.12)') 'GENESIS_NONBOND_VIRIAL: CPU reference pass ',ratio
    active_mode=0
  end subroutine
end module
