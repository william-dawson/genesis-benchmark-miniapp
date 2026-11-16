#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module hgr_migration_mod
  use iso_c_binding,only:c_int,c_int16_t,c_double
  use constants_mod,only:wp,wip
  use cudafor
  use sp_domain_str_mod
  use sp_boundary_str_mod
  use sp_constraints_str_mod
  use mpi_parallel_mod
  implicit none
  private
  public::hgr_latch_c,hgr_outgoing,hgr_check,hgr_verify_c
  interface
    integer(c_int) function enabled_c() bind(C,name='genesis_hgr_enabled')
      import c_int
    end function
    integer(c_int) function hgr_verify_c() bind(C,name='genesis_hgr_verify')
      import c_int
    end function
    subroutine hgr_latch_c(x,v,groups,counts,dist,width,cells,group_width,sites,connect) bind(C,name='genesis_hgr_latch')
      import c_int,wip
      real(wip),device,intent(in)::x(*),v(*),dist(*)
      integer(c_int),device,intent(in)::groups(*),counts(*)
      integer(c_int),value::width,cells,group_width,sites,connect
    end subroutine
    integer(c_int) function outgoing_c(counts,groups,ids,solute,cls,charge,mass,map,box,origin,cell_size,grid,atoms, &
        stay,move,stay_real,move_real,stay_integer,move_integer,stay_capacity,move_capacity) bind(C,name='genesis_hgr_outgoing')
      import c_int,c_int16_t,c_double,wp,wip
      integer(c_int),intent(in)::counts(*),groups(*),ids(*),solute(*),cls(*),grid(*)
      integer(c_int16_t),intent(in)::map(*)
      real(wp),intent(in)::charge(*)
      real(wip),intent(in)::mass(*)
      real(c_double),intent(in)::box(*),origin(*),cell_size(*)
      integer(c_int),value::atoms,stay_capacity,move_capacity
      integer(c_int),intent(out)::stay(*),move(*),stay_integer(*),move_integer(*)
      real(wip),intent(out)::stay_real(*),move_real(*)
    end function
    subroutine check_c(stay,move,stay_real,move_real,stay_integer,move_integer,stay_capacity,move_capacity) &
        bind(C,name='genesis_hgr_check')
      import c_int,wip
      integer(c_int),intent(in)::stay(*),move(*),stay_integer(*),move_integer(*)
      real(wip),intent(in)::stay_real(*),move_real(*)
      integer(c_int),value::stay_capacity,move_capacity
    end subroutine
  end interface
contains
  subroutine scatter(stay_real,move_real,stay_integer,move_integer,stay,move,packed,integer,ranges, &
      cells,connect,sites,stay_capacity,move_capacity) bind(C,name='genesis_hgr_scatter')
    integer(c_int),value::cells,connect,sites,stay_capacity,move_capacity
    real(wip),intent(out)::stay_real(9*sites,stay_capacity,connect,cells),move_real(9*sites,move_capacity,connect,cells)
    integer(c_int),intent(out)::stay_integer(3*sites,stay_capacity,connect,cells)
    integer(c_int),intent(out)::move_integer(3*sites,move_capacity,connect,cells),stay(connect,cells),move(connect,cells)
    real(wip),intent(in)::packed(9*sites,*)
    integer(c_int),intent(in)::integer(3*sites,*),ranges(2,2*connect*cells)
    integer::cell,j,b,first,last,count
    !$omp parallel do private(j,b,first,last,count) schedule(static)
    do cell=1,cells
      do j=1,connect
        b=j+connect*(cell-1)
        first=ranges(1,b)+1;last=ranges(2,b);count=0
        if(first>0)count=last-first+1
        stay(j,cell)=count
        if(count>0)then
          stay_real(:,1:count,j,cell)=packed(:,first:last)
          stay_integer(:,1:count,j,cell)=integer(:,first:last)
        end if
        b=b+connect*cells
        first=ranges(1,b)+1;last=ranges(2,b);count=0
        if(first>0)count=last-first+1
        move(j,cell)=count
        if(count>0)then
          move_real(:,1:count,j,cell)=packed(:,first:last)
          move_integer(:,1:count,j,cell)=integer(:,first:last)
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine
  logical function hgr_outgoing(boundary,constraints,domain)
    type(s_boundary),intent(in)::boundary
    type(s_constraints),intent(inout)::constraints
    type(s_domain),intent(inout)::domain
    real(8)::box(3),origin(3),cell_size(3)
    integer::grid(3)
    hgr_outgoing=.false.
    if(enabled_c()==0)return
    if(nproc_country/=1.or.domain%num_cell_boundary/=0) error stop 'GPU HGr requires one rank without halos'
    if(size(domain%id_l2g,1)/=size(domain%coord,2))error stop 'HGr ID pitch mismatch'
    box=[boundary%box_size_x,boundary%box_size_y,boundary%box_size_z]
    origin=[boundary%origin_x,boundary%origin_y,boundary%origin_z]
    cell_size=[boundary%cell_size_x,boundary%cell_size_y,boundary%cell_size_z]
    grid=[boundary%num_cells_x,boundary%num_cells_y,boundary%num_cells_z]
    hgr_outgoing=outgoing_c(constraints%HGr_local(1,1),constraints%HGr_bond_list(1,1,1,1),domain%id_l2g(1,1), &
      domain%id_l2g_solute(1,1),domain%atom_cls_no(1,1),domain%charge(1,1),domain%mass(1,1),domain%cell_g2l, &
      box,origin,cell_size,grid,int(domain%num_atom_all),constraints%HGr_stay(1,1),constraints%HGr_move(1,1), &
      constraints%HGr_stay_real(1,1,1,1),constraints%HGr_move_real(1,1,1,1), &
      constraints%HGr_stay_int(1,1,1,1),constraints%HGr_move_int(1,1,1,1), &
      size(constraints%HGr_stay_int,2),size(constraints%HGr_move_int,2))==1
  end function
  subroutine hgr_check(constraints)
    type(s_constraints),intent(in)::constraints
    call check_c(constraints%HGr_stay(1,1),constraints%HGr_move(1,1), &
      constraints%HGr_stay_real(1,1,1,1),constraints%HGr_move_real(1,1,1,1), &
      constraints%HGr_stay_int(1,1,1,1),constraints%HGr_move_int(1,1,1,1), &
      size(constraints%HGr_stay_int,2),size(constraints%HGr_move_int,2))
  end subroutine
end module
