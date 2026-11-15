#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module remap_coordinate_mod
  use iso_c_binding,only:c_int,c_double
  use cudafor
  use sp_domain_str_mod
  implicit none
  private
  public::remap_mode,remap_allowed,remap_publish,remap_record_c,remap_consume_c
  interface
    integer(c_int) function remap_mode() bind(C,name='genesis_remap_mode')
      import c_int
    end function
    integer(c_int) function allowed_c() bind(C,name='genesis_remap_allowed')
      import c_int
    end function
    subroutine prepare_c(coord,velocity,gpu_velocity,translation,charge,counts,starts) &
        bind(C,name='genesis_coord_prepare_remap')
      import c_double,c_int
      real(c_double),intent(in)::coord(*),velocity(*),translation(*),charge(*)
      real(c_double),device,intent(in)::gpu_velocity(*)
      integer(c_int),intent(in)::counts(*),starts(*)
    end subroutine
    subroutine remap_record_c(consumer) bind(C,name='genesis_remap_record')
      import c_int
      integer(c_int),value::consumer
    end subroutine
    subroutine remap_consume_c() bind(C,name='genesis_remap_consume')
    end subroutine
  end interface
contains
  logical function remap_allowed()
    remap_allowed=allowed_c()/=0
  end function
  subroutine remap_publish(domain,velocity)
    type(s_domain),intent(in)::domain
    real(8),device,intent(in)::velocity(:,:,:)
    call prepare_c(domain%coord(1,1,1),domain%velocity(1,1,1),velocity,domain%trans_vec(1,1,1), &
      domain%charge(1,1),domain%num_atom,domain%start_atom)
  end subroutine
end module
