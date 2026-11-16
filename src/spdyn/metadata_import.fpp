#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module metadata_import_mod
  use iso_c_binding,only:c_double,c_int
  use constants_mod,only:wip
  use cudafor
  use sp_domain_str_mod
  use sp_constraints_str_mod
  implicit none
  private
  public::metadata_import
  interface
    integer(c_int) function load_c(m,im,na,nw,wl,ns,hc,hl,hd,hm,him,hna,hnw,hwl,hns,hhc,hhl,hhd, &
        width,cells,wh,wc,connect,gh,gc,dh,dc,eligible) bind(C,name='genesis_metadata_import')
      import wip,c_int
      real(wip),device::m(*),im(*),hd(*)
      integer(c_int),device::na(*),nw(*),wl(*),ns(*),hc(*),hl(*)
      real(wip),intent(in)::hm(*),him(*),hhd(*)
      integer(c_int),intent(in)::hna(*),hnw(*),hwl(*),hns(*),hhc(*),hhl(*)
      integer(c_int),value::width,cells,wh,wc,connect,gh,gc,dh,dc,eligible
    end function
  end interface
contains
  logical function metadata_import(domain,constraints,m,im,na,nw,wl,ns,hc,hl,hd,eligible) result(loaded)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(wip),device,intent(inout)::m(:,:),im(:,:),hd(:,:,:,:)
    integer,device,intent(inout)::na(:),nw(:),wl(:,:,:),ns(:),hc(:,:),hl(:,:,:,:)
    logical,intent(in)::eligible
    loaded=load_c(m,im,na,nw,wl,ns,hc,hl,hd,domain%mass,domain%inv_mass,domain%num_atom, &
      domain%num_water,domain%water_list,constraints%No_HGr,constraints%HGr_local, &
      constraints%HGr_bond_list,constraints%HGr_bond_dist,size(m,1),size(m,2),size(wl,1),size(wl,2), &
      size(hc,1),size(hl,1),size(hl,2),size(hd,1),size(hd,2),merge(1,0,eligible))/=0
  end function
end module
