!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_update_domain_mod
!> @brief   library subroutine used for integrator
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_update_domain_mod
#ifdef USE_GPU
  use translation_rebuild_mod,only:translation_rebuild_finish
#endif

  use sp_domain_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_communicate_mod
  use sp_migration_mod
  use sp_dynamics_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use maintenance_probe_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: domain_interaction_update

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_interaction_update
  !> @brief        update interactions with MD
  !! @authors      JJ
  !! @param[in]    istep         : step number
  !! @param[in]    update_period : update period for pairlist
  !! @param[inout] domain        : domain information
  !! @param[inout] enefunc       : potential energy functions information
  !! @param[inout] pairlist      : pairlist information
  !! @param[inout] boundary      : boundary condition information
  !! @param[inout] constraints   : constraints information
  !! @param[inout] comm          : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_interaction_update(istep, update_period, domain, enefunc, &
                                       pairlist, boundary, constraints, comm)

    ! formal arguments
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: update_period
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: ncell, nb, water_atom


    call maintenance_probe_begin(17)
    ncell = domain%num_cell_local
    nb    = domain%num_cell_boundary
    pairlist%univ_update = 0

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    if (mod(istep,update_period) == 0) then

#ifdef DEBUG
      if (main_rank) then
        write(MsgOut,'(a,3f15.8)') &
             "Debuging > bsize_[x,y,z]", &
             boundary%box_size_x,boundary%box_size_y,boundary%box_size_z
        write(MsgOut,'(a,2i8)')    &
             "Debuging > ncell, nb",ncell,nb
      end if
#endif

      pairlist%univ_update = 1

      if (istep == update_period) then

    call maintenance_probe_begin(2)
        call alloc_domain(domain, DomainPtlMove, ncell+nb, ncell, 1)
    call maintenance_probe_end(2)
    call maintenance_probe_begin(2)
        call alloc_domain(domain, DomainWaterMove, ncell+nb, ncell, &
                          water_atom)
    call maintenance_probe_end(2)

    call maintenance_probe_begin(3)
        call alloc_constraints(constraints, ConstraintsHGroupMove, &
                               constraints%connect, ncell+nb)
    call maintenance_probe_end(3)

    call maintenance_probe_begin(4)
        call alloc_enefunc(enefunc, EneFuncBondCell, ncell, ncell+nb)
    call maintenance_probe_end(4)

      end if

      call timer(TimerUpdate, TimerOn)
    call maintenance_probe_begin(5)
      call update_outgoing_atom (boundary, constraints, domain)
    call maintenance_probe_end(5)
    call maintenance_probe_begin(6)
      call update_outgoing_HGr  (boundary, constraints, domain)
    call maintenance_probe_end(6)
    call maintenance_probe_begin(7)
      call update_outgoing_water(water_atom, boundary, domain)
    call maintenance_probe_end(7)

      call timer(TimerComm3, TimerOn)
    call maintenance_probe_begin(8)
      call communicate_constraints(domain, comm, constraints)
    call maintenance_probe_end(8)
      call timer(TimerComm3, TimerOff)

    call maintenance_probe_begin(9)
      call update_incoming_atom (constraints, domain)
    call maintenance_probe_end(9)
    call maintenance_probe_begin(10)
      call update_incoming_HGr  (constraints, domain)
    call maintenance_probe_end(10)
    call maintenance_probe_begin(11)
      call update_incoming_water(water_atom, domain)
    call maintenance_probe_end(11)

    call maintenance_probe_begin(12)
      call update_cell_size_constraints    (domain, comm, constraints)
    call maintenance_probe_end(12)
    call maintenance_probe_begin(13)
      call update_cell_boundary(domain, comm, boundary, constraints)
#ifdef USE_GPU
      call translation_rebuild_finish(domain,comm,boundary,constraints)
#endif
    call maintenance_probe_end(13)

    call maintenance_probe_begin(14)
      call update_enefunc(domain, comm, enefunc, constraints)
    call maintenance_probe_end(14)

      call timer(TimerUpdate, TimerOff)

      call timer(TimerPairList, TimerOn)

      if (enefunc%pairlist_check) then
    call maintenance_probe_begin(15)
        call update_pairlist_pbc_check(enefunc, domain, pairlist)
    call maintenance_probe_end(15)
      else
    call maintenance_probe_begin(16)
        call update_pairlist_pbc(enefunc, domain, pairlist)
    call maintenance_probe_end(16)
      end if

      call timer(TimerPairList, TimerOff)

    end if

    call maintenance_probe_end(17)
    return

  end subroutine domain_interaction_update
 
end module sp_update_domain_mod
