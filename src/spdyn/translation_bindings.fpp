#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module translation_residency_mod
  use iso_c_binding,only:c_int,c_double
  use constants_mod
  use sp_domain_str_mod
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  implicit none
  private
  public::translation_mode,translation_reference,translation_scaled,translation_configure,translation_refresh
  real(wp),allocatable,save::reference(:,:,:)
  interface
    integer(c_int) function translation_mode() bind(C,name='genesis_translation_mode')
      import c_int
    end function
    subroutine publish_c(values) bind(C,name='genesis_translation_publish')
      import wp
      real(wp),intent(in)::values(*)
    end subroutine
    subroutine refresh_c(values) bind(C,name='genesis_translation_refresh')
      import wp
      real(wp),intent(inout)::values(*)
    end subroutine
  end interface
contains
  subroutine translation_configure(domain)
    type(s_domain),intent(in)::domain
    if(translation_mode()/=2) return
    if(.not.allocated(reference)) allocate(reference,mold=domain%trans_vec)
    if(any(shape(reference)/=shape(domain%trans_vec))) error stop 'Translation reference shape changed'
    reference=domain%trans_vec
  end subroutine

  subroutine translation_reference(domain,scale)
    type(s_domain),intent(in)::domain
    real(wip),intent(in)::scale(3)
    integer::cell,atom
    if(translation_mode()/=2) return
    if(.not.allocated(reference)) error stop 'Missing translation reference'
    !$omp parallel do private(atom) schedule(static)
    do cell=1,domain%num_cell_local
      do atom=1,domain%num_atom(cell)
        reference(:,atom,cell)=reference(:,atom,cell)*scale
      end do
    end do
    !$omp end parallel do
  end subroutine

  subroutine translation_scaled(domain)
    type(s_domain),intent(inout)::domain
    if(translation_mode()==0) return
    if(translation_mode()==2) then
      call publish_c(reference(1,1,1))
      domain%trans_vec=ieee_value(0.0_wp,ieee_quiet_nan)
      write(*,'(a)') 'GENESIS_TRANSLATION_VERIFY: CPU mirror poisoned'
    else
      call publish_c(domain%trans_vec(1,1,1))
    end if
  end subroutine

  subroutine translation_refresh(domain)
    type(s_domain),intent(inout)::domain
    if(translation_mode()/=0) call refresh_c(domain%trans_vec(1,1,1))
  end subroutine
end module
