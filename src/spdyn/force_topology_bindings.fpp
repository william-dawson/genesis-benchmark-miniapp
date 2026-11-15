module force_topology_bindings_mod
 use iso_c_binding
 use cudafor
 implicit none
 integer,save::host_live_bond=0,host_live_angle=0,host_live_dihed=0,host_live_improper=0,host_live_nb14=0,host_live_corr=0
 interface
  subroutine topology_fields(kind,cells,expected,capacity,offset,owner,local,counts,lists,pbc,period,r0,r1,r2,r3) &
      bind(C,name='genesis_force_topology_fields')
   import c_int,c_double
   integer(c_int),value::kind,cells,expected,capacity
   integer(c_int),device::offset(*),owner(*),local(*),counts(*),lists(*),pbc(*),period(*)
   real(c_double),device::r0(*),r1(*),r2(*),r3(*)
  end subroutine
  subroutine topology_exclusion(kind,cells,expected,capacity,offset,owner,local,counts,starts,lists,move,host_move,pitch) &
      bind(C,name='genesis_force_topology_exclusion')
   import c_int,c_int8_t
   integer(c_int),value::kind,cells,expected,capacity,pitch
   integer(c_int),device::offset(*),owner(*),local(*),counts(*),starts(*),lists(*)
   integer(c_int8_t),device::move(*)
   integer(c_int8_t),intent(in)::host_move(*)
  end subroutine
  subroutine topology_start(stream) bind(C,name='genesis_force_topology_start')
   import cuda_stream_kind
   integer(kind=cuda_stream_kind),value::stream
  end subroutine
  integer(c_int) function topology_bypass() bind(C,name='genesis_force_topology_bypass')
   import c_int
  end function
  subroutine topology_check(stream) bind(C,name='genesis_force_topology_check')
   import cuda_stream_kind
   integer(kind=cuda_stream_kind),value::stream
  end subroutine
 end interface
end module
