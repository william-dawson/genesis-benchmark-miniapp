module maintenance_probe_mod
  implicit none
  private
  public::maintenance_probe_begin,maintenance_probe_end,maintenance_probe_reset,maintenance_probe_dump
  integer,parameter::stages=33
  character(64),parameter::names(stages)=[character(64):: &
      'production_md', &
      'alloc_domain', &
      'alloc_constraints', &
      'alloc_enefunc', &
      'update_outgoing_atom', &
      'update_outgoing_hgr', &
      'update_outgoing_water', &
      'communicate_constraints', &
      'update_incoming_atom', &
      'update_incoming_hgr', &
      'update_incoming_water', &
      'update_cell_size_constraints', &
      'update_cell_boundary', &
      'update_enefunc', &
      'update_pairlist_pbc_check', &
      'update_pairlist_pbc', &
      'domain_driver', &
      'update_outgoing_enefunc_bond', &
      'update_outgoing_enefunc_angl', &
      'update_outgoing_enefunc_dihe', &
      'update_outgoing_enefunc_rb_dihe', &
      'update_outgoing_enefunc_impr', &
      'update_outgoing_enefunc_restraint', &
      'update_outgoing_enefunc_fitting', &
      'communicate_bond', &
      'update_incoming_enefunc_bond', &
      'update_incoming_enefunc_angl', &
      'update_incoming_enefunc_dihe', &
      'update_incoming_enefunc_rb_dihe', &
      'update_incoming_enefunc_impr', &
      'update_incoming_enefunc_restraint', &
      'update_incoming_enefunc_fitting', &
      'count_nonb_excl']
  logical,save::collect=.false.
  integer,save::depth=0,stack(32),calls(stages)=0
  real(8),save::started(32),children(32),wall(stages)=0,self(stages)=0
contains
  real(8) function now()
    integer(8)::ticks,rate
    call system_clock(ticks,rate)
    now=real(ticks,8)/real(rate,8)
  end function
  subroutine maintenance_probe_reset()
    integer::status
    character(16)::value
    call get_environment_variable('GENESIS_MAINT_PROBE',value,status=status)
    collect=status==0.and.trim(value)=='1'
    depth=0;calls=0;wall=0;self=0
  end subroutine
  subroutine maintenance_probe_begin(stage)
    integer,intent(in)::stage
    if(.not.collect) return
    depth=depth+1
    if(depth>32) error stop 'Maintenance probe nesting overflow'
    stack(depth)=stage;started(depth)=now();children(depth)=0
  end subroutine
  subroutine maintenance_probe_end(stage)
    integer,intent(in)::stage
    real(8)::elapsed
    if(.not.collect) return
    if(depth<1) error stop 'Maintenance probe nesting underflow'
    if(stack(depth)/=stage) error stop 'Maintenance probe nesting mismatch'
    elapsed=now()-started(depth)
    calls(stage)=calls(stage)+1;wall(stage)=wall(stage)+elapsed
    self(stage)=self(stage)+elapsed-children(depth)
    depth=depth-1
    if(depth>0) children(depth)=children(depth)+elapsed
  end subroutine
  subroutine maintenance_probe_dump()
    integer::unit,stage
    if(.not.collect) return
    if(depth/=0) error stop 'Unclosed maintenance probe region'
    collect=.false.
    open(newunit=unit,file='maintenance.csv',status='replace',action='write')
    write(unit,'(a)') 'stage,calls,inclusive_s,exclusive_s'
    do stage=1,stages
      if(calls(stage)==0) cycle
      write(unit,'(a,a,i0,2(a,es24.16))') trim(names(stage)),',',calls(stage),',',wall(stage),',',self(stage)
    end do
    close(unit)
  end subroutine
end module
