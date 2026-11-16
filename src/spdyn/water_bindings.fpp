#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module water_migration_mod
  use iso_c_binding,only:c_int,c_int16_t,c_double
  use constants_mod,only:wip
  use cudafor
  use sp_domain_str_mod
  use sp_boundary_str_mod
  use mpi_parallel_mod
  implicit none
  private
  public::water_latch_c,water_outgoing,water_check,water_verify_c
  interface
    integer(c_int) function water_enabled_c() bind(C,name='genesis_water_enabled')
      import c_int
    end function
    integer(c_int) function water_verify_c() bind(C,name='genesis_water_verify')
      import c_int
    end function
    subroutine water_latch_c(x,v,water,counts,width,cells,water_width) bind(C,name='genesis_water_latch')
      import c_int,wip
      real(wip),device,intent(in)::x(*),v(*)
      integer(c_int),device,intent(in)::water(*),counts(*)
      integer(c_int),value::width,cells,water_width
    end subroutine
    integer(c_int) function outgoing_c(counts,ids,map,box,origin,cell_size,grid,stay,move, &
        stay_real,move_real,stay_ids,move_ids,stay_capacity,move_capacity) bind(C,name='genesis_water_outgoing')
      import c_int,c_int16_t,c_double,wip
      integer(c_int),intent(in)::counts(*),ids(*),grid(*)
      integer(c_int16_t),intent(in)::map(*)
      real(c_double),intent(in)::box(*),origin(*),cell_size(*)
      integer(c_int),intent(out)::stay(*),move(*),stay_ids(*),move_ids(*)
      real(wip),intent(out)::stay_real(*),move_real(*)
      integer(c_int),value::stay_capacity,move_capacity
    end function
    subroutine check_c(stay,move,stay_real,move_real,stay_ids,move_ids,stay_capacity,move_capacity) &
        bind(C,name='genesis_water_check')
      import c_int,wip
      integer(c_int),intent(in)::stay(*),move(*),stay_ids(*),move_ids(*)
      real(wip),intent(in)::stay_real(*),move_real(*)
      integer(c_int),value::stay_capacity,move_capacity
    end subroutine
  end interface
contains
  subroutine scatter(stay_real,move_real,stay_ids,move_ids,stay,move,packed,packed_ids,ranges, &
      cells,stay_capacity,move_capacity) bind(C,name='genesis_water_scatter')
    integer(c_int),value::cells,stay_capacity,move_capacity
    real(wip),intent(out)::stay_real(18,stay_capacity,cells),move_real(18,move_capacity,cells)
    integer(c_int),intent(out)::stay_ids(3,stay_capacity,cells),move_ids(3,move_capacity,cells),stay(cells),move(cells)
    real(wip),intent(in)::packed(18,*)
    integer(c_int),intent(in)::packed_ids(3,*),ranges(2,2*cells)
    integer::cell,first,last,count
    !$omp parallel do private(first,last,count) schedule(static)
    do cell=1,cells
      first=ranges(1,cell)+1;last=ranges(2,cell)
      count=0
      if(first>0)count=last-first+1
      stay(cell)=count
      if(count>0)then
        stay_real(:,1:count,cell)=packed(:,first:last)
        stay_ids(:,1:count,cell)=packed_ids(:,first:last)
      end if
      first=ranges(1,cells+cell)+1;last=ranges(2,cells+cell)
      count=0
      if(first>0)count=last-first+1
      move(cell)=count
      if(count>0)then
        move_real(:,1:count,cell)=packed(:,first:last)
        move_ids(:,1:count,cell)=packed_ids(:,first:last)
      end if
    end do
    !$omp end parallel do
  end subroutine
  logical function water_outgoing(water_atom,boundary,domain)
    integer,intent(in)::water_atom
    type(s_boundary),intent(in)::boundary
    type(s_domain),intent(inout)::domain
    real(8)::box(3),origin(3),cell_size(3)
    integer::grid(3)
    water_outgoing=.false.
    if(water_enabled_c()==0)return
    if(nproc_country/=1.or.domain%num_cell_boundary/=0.or.water_atom/=3) &
      error stop 'GPU water migration requires one rank, no halos and three-site waters'
    if(size(domain%id_l2g,1)/=size(domain%coord,2)) error stop 'GPU water ID pitch mismatch'
    box=[boundary%box_size_x,boundary%box_size_y,boundary%box_size_z]
    origin=[boundary%origin_x,boundary%origin_y,boundary%origin_z]
    cell_size=[boundary%cell_size_x,boundary%cell_size_y,boundary%cell_size_z]
    grid=[boundary%num_cells_x,boundary%num_cells_y,boundary%num_cells_z]
    water_outgoing=outgoing_c(domain%num_water,domain%id_l2g(1,1),domain%cell_g2l,box,origin,cell_size,grid, &
      domain%water%stay,domain%water%move,domain%water%stay_real(1,1,1),domain%water%move_real(1,1,1), &
      domain%water%stay_integer(1,1,1),domain%water%move_integer(1,1,1), &
      size(domain%water%stay_integer,2),size(domain%water%move_integer,2))==1
  end function
  subroutine water_check(domain)
    type(s_domain),intent(in)::domain
    call check_c(domain%water%stay,domain%water%move,domain%water%stay_real(1,1,1), &
      domain%water%move_real(1,1,1),domain%water%stay_integer(1,1,1),domain%water%move_integer(1,1,1), &
      size(domain%water%stay_integer,2),size(domain%water%move_integer,2))
  end subroutine
end module
