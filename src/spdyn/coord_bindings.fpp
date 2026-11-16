#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module coordinate_pack_mod
  use iso_c_binding,only:c_int,c_double,c_intptr_t
  use constants_mod,only:wp,wip
  use cudafor
  use sp_domain_str_mod
  use translation_residency_mod,only:translation_configure,translation_refresh
  use mpi_parallel_mod
  implicit none
  private
  public::coordinate_enabled,coordinate_active,coordinate_resident,coordinate_verify
  public::coordinate_configure,coordinate_check,coordinate_cpu_required
  public::coordinate_pack_c,coordinate_invalidate_c,coordinate_copy_c
  public::coordinate_pme,coordinate_host_refreshed_c,coordinate_require_host_c
  interface
    integer(c_int) function enabled_c() bind(C,name='genesis_coord_enabled')
      import c_int
    end function
    integer(c_int) function active_c() bind(C,name='genesis_coord_active')
      import c_int
    end function
    integer(c_int) function resident_c() bind(C,name='genesis_coord_resident')
      import c_int
    end function
    integer(c_int) function verify_c() bind(C,name='genesis_coord_verify')
      import c_int
    end function
    subroutine coordinate_invalidate_c() bind(C,name='genesis_coord_invalidate')
    end subroutine
    integer(c_int) function pme_c() bind(C,name='genesis_coord_pme')
      import c_int
    end function
    subroutine coordinate_host_refreshed_c() bind(C,name='genesis_coord_host_refreshed')
    end subroutine
    subroutine coordinate_require_host_c() bind(C,name='genesis_coord_require_host')
    end subroutine
    subroutine configure_c(values,translation,charge,counts,starts,width,cells,atoms) bind(C,name='genesis_coord_configure')
      import wip,c_int,wp
      real(wip),device,intent(in)::values(*)
      real(wp),intent(in)::translation(*),charge(*)
      integer(c_int),intent(in)::counts(*),starts(*)
      integer(c_int),value::width,cells,atoms
    end subroutine
    subroutine coordinate_pack_c(x,y,z,short_step) bind(C,name='genesis_coord_pack')
      import c_double,c_int
      real(c_double),value::x,y,z
      integer(c_int),value::short_step
    end subroutine
    integer(c_int) function coordinate_copy_c(values,atoms,components,stream,consumer) bind(C,name='genesis_coord_copy')
      import wp,c_int,c_intptr_t
      real(wp),device::values(*)
      integer(c_int),value::atoms,components,consumer
      integer(c_intptr_t),value::stream
    end function
    integer(c_int) function cpu_required_c(virial,cpu_calc) bind(C,name='genesis_coord_cpu_required')
      import c_int
      integer(c_int),value::virial,cpu_calc
    end function
    subroutine check_c(values,translation,charge,counts,starts) bind(C,name='genesis_coord_check')
      import wip,wp,c_int
      real(wip),intent(in)::values(*)
      real(wp),intent(in)::translation(*),charge(*)
      integer(c_int),intent(in)::counts(*),starts(*)
    end subroutine
  end interface
contains
  logical function coordinate_enabled()
    coordinate_enabled=enabled_c()==1
  end function
  logical function coordinate_active()
    coordinate_active=active_c()==1
  end function
  logical function coordinate_resident()
    coordinate_resident=resident_c()==1
  end function
  integer function coordinate_verify()
    coordinate_verify=verify_c()
  end function
  logical function coordinate_pme()
    coordinate_pme=pme_c()==1
  end function
  logical function coordinate_cpu_required(virial,cpu_calc)
    logical,intent(in)::virial,cpu_calc
    coordinate_cpu_required=cpu_required_c(merge(1,0,virial),merge(1,0,cpu_calc))==1
  end function
  subroutine coordinate_configure(domain,values)
    type(s_domain),intent(in)::domain
    real(wip),device,intent(in)::values(:,:,:)
    if(.not.coordinate_enabled()) return
    if(nproc_country/=1.or.domain%num_cell_boundary/=0) error stop 'Resident coordinates require one rank without halos'
    call configure_c(values,domain%trans_vec(1,1,1),domain%charge(1,1),domain%num_atom,domain%start_atom, &
      size(values,2),size(values,3),domain%num_atom_domain)
    call translation_configure(domain)
  end subroutine
  subroutine coordinate_check(domain)
    type(s_domain),intent(inout)::domain
    if(coordinate_verify()/=0) call translation_refresh(domain)
    call check_c(domain%coord(1,1,1),domain%trans_vec(1,1,1),domain%charge(1,1),domain%num_atom,domain%start_atom)
  end subroutine
end module
