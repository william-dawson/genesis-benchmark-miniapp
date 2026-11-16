#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS	:= ../lib .

include $(top_srcdir)/src/mk/common-libonly.am

bin_PROGRAMS = spdyn

spdyn_SOURCES = \
	domain_masks.fpp \
	exclusion_reset.fpp \
	exclusion_oracle.fpp \
	sp_domain_str.fpp \
	sp_energy_str.fpp \
	sp_enefunc_str.fpp \
	sp_pairlist_str.fpp \
	sp_boundary_str.fpp \
	sp_constraints_str.fpp \
	sp_restraints_str.fpp \
	sp_ensemble_str.fpp \
	sp_dynvars_str.fpp \
	sp_dynamics_str.fpp \
	sp_output_str.fpp \
	incoming_oracle.fpp \
	sp_migration.fpp \
	sp_communicate.fpp \
	sp_pairlist_generic.fpp \
	sp_pairlist_fugaku.fpp \
	pair_sparse.fpp \
	pair_host.fpp \
	pair_sort.fpp \
	pair_pack.c \
	sp_pairlist_gpu.fpp \
	sp_pairlist.fpp \
	sp_boundary.fpp \
	sp_energy_efield.fpp \
	sp_energy_bonds.fpp \
	sp_energy_angles.fpp \
	sp_energy_dihedrals.fpp \
	sp_energy_nonbond_generic.fpp \
	sp_energy_nonbond_fugaku.fpp \
	sp_energy_nonbond_intel.fpp \
	sp_energy_nonbond_gpu.fpp \
	sp_energy_table_linear.fpp \
	sp_energy_14_notable.fpp \
	sp_energy_nonbond_pme_table.fpp \
	sp_energy_nonbond_pme_notable.fpp \
	sp_energy_nonbond_ljpme.fpp \
	sp_energy_table_linear_bondcorr.fpp \
	sp_energy_table_cubic.fpp \
	sp_energy_pme_opt_1dalltoall.fpp \
	sp_energy_pme_noopt_1dalltoall.fpp \
	sp_energy_pme_opt_2dalltoall.fpp \
	sp_energy_pme_noopt_2dalltoall.fpp \
	sp_energy_pme_opt_slab.fpp \
	sp_energy_pme_noopt_slab.fpp \
	sp_energy_pme.fpp \
	sp_energy_nonbonds.fpp \
	sp_enefunc_fit.fpp \
	sp_energy_restraints.fpp \
	sp_energy.fpp \
	sp_enefunc_table.fpp \
	sp_enefunc_restraints.fpp \
	sp_enefunc_localres.fpp \
	sp_enefunc_charmm.fpp \
	sp_enefunc_amber.fpp \
	sp_enefunc_gromacs.fpp \
	sp_enefunc.fpp \
	sp_constraints.fpp \
	sp_restraints.fpp \
	sp_ensemble.fpp \
	sp_domain.fpp \
	sp_dynvars.fpp \
	sp_assign_velocity.fpp \
	maintenance_probe.fpp \
	sp_update_domain.fpp \
	sp_input.fpp \
	sp_output.fpp \
	sp_md_leapfrog.fpp \
	sp_md_vverlet.fpp \
	sp_md_respa.fpp \
	sp_dynamics.fpp \
	sp_restart.fpp \
	sp_control.fpp \
	sp_setup_spdyn.fpp \
	sp_setup_mpi.fpp 

if CUDAGPU
PME_CUFFT_LIBS = -lcufft
STEP_CUDA_LIBS = -cuda
gpu_thermo.o: FCFLAGS += -cuda -gpu=cc100,keep,ptxinfo,nofma
gpu_thermo.o: sp_constraints.o sp_domain_str.o sp_constraints_str.o precision_rounding.o
precision_rounding.o: sp_domain_str.o sp_constraints_str.o
water_bindings.o: FCFLAGS += -cuda
water_bindings.o: sp_domain_str.o sp_boundary_str.o
gpu_step.o sp_migration.o: water_bindings.o
hgr_bindings.o: FCFLAGS += -cuda
hgr_bindings.o: sp_domain_str.o sp_boundary_str.o sp_constraints_str.o
gpu_step.o sp_migration.o: hgr_bindings.o
coord_bindings.o: FCFLAGS += -cuda
remap_coordinate.o: FCFLAGS += -cuda
remap_coordinate.o: sp_domain_str.o
final_reductions.o: FCFLAGS += -cuda
gpu_step.o gpu_thermo.o: final_reductions.o
incoming_gpu.o: incoming_oracle.o sp_domain_str.o sp_constraints_str.o sp_enefunc_str.o
sp_migration.o sp_communicate.o translation_rebuild.o: incoming_gpu.o
bonded_incoming.o: bonded_migration.o sp_domain_str.o sp_enefunc_str.o
sp_migration.o: bonded_incoming.o
gpu_step.o: remap_coordinate.o
nonbond_virial.o: sp_domain_str.o force_assembly.o
sp_energy.o sp_energy_nonbond_gpu.o: nonbond_virial.o
maintenance_state.o: FCFLAGS += -cuda
maintenance_state.o: sp_domain_str.o
gpu_step.o: maintenance_state.o
translation_rebuild.o: translation_bindings.o sp_communicate.o sp_domain_str.o sp_constraints_str.o sp_boundary_str.o
sp_update_domain.o sp_md_respa.o: translation_rebuild.o
translation_bindings.o: sp_domain_str.o
coord_bindings.o: sp_domain_str.o translation_bindings.o
gpu_step.o: coord_bindings.o
bonded_force.o sp_energy_nonbond_gpu.o sp_md_respa.o: coord_bindings.o
gpu_step.o: gpu_thermo.o
gpu_step.o: FCFLAGS += -cuda -gpu=cc100,keep,ptxinfo,nofma
gpu_step.o: sp_domain_str.o sp_constraints_str.o sp_dynamics_str.o sp_ensemble_str.o sp_dynvars_str.o sp_boundary_str.o sp_enefunc_str.o
sp_md_respa.o: gpu_step.o
gpu_step.o: force_assembly.o sp_assign_velocity.o
force_assembly.o: sp_domain_str.o
bonded_force.o: FCFLAGS += -cuda -gpu=cc100,keep,ptxinfo,nofma
bonded_force.o: gpu_step.o sp_enefunc_str.o sp_energy_str.o sp_pairlist_str.o force_assembly.o
sp_energy.o: force_assembly.o bonded_force.o
spdyn_SOURCES += \
	bonded_force.fpp \
	force_assembly.fpp \
	force_totals_gpu.cu \
	force_topology_bindings.fpp \
	metadata_import.fpp \
	state_cycle.fpp \
	state_cycle_gpu.cu \
	observables.fpp \
	observables_gpu.cu \
	metadata_import_gpu.cu \
	force_topology_gpu.cu \
	precision_types.h \
	migration_coordinates.h \
	migration_views.h \
	bonded_incoming_view.h \
	long_force_gpu.cu \
	gpu_state_ops.cu \
	gpu_force_assembly.cu \
	remap_coordinate.fpp \
	remap_coordinate_gpu.cu \
	nonbond_mirror_gpu.cu \
	nonbond_virial.fpp \
	nonbond_virial_gpu.cu \
	final_reductions.fpp \
	final_reductions_gpu.cu \
	incoming_gpu.fpp \
	incoming_gpu_ops.cu \
	atom_outgoing_gpu.cu \
	migration_resident_gpu.cu \
	bonded_incoming.fpp \
	bonded_incoming_gpu.cu \
	bonded_incoming_hbm.cu \
	bonded_incoming_cooperative.cu \
	improper_outgoing.fpp \
	improper_outgoing_gpu.cu \
	maintenance_state.fpp \
	maintenance_state_gpu.cu \
	translation_rebuild.fpp \
	translation_bindings.fpp \
	translation_gpu.cu \
	coord_bindings.fpp \
	gpu_coord_pack.cu \
	hgr_bindings.fpp \
	hgr_migration.cu \
	water_bindings.fpp \
	water_migration.cu \
	precision_rounding.fpp \
	gpu_thermo.fpp \
	gpu_step.fpp \
	gpu_sp_energy.cu \
	nonbond_inter_staged.cu \
	nonbond_inter_staged.h \
	gpu_pme.cu \
	pair_sparse_gpu.cu \
	pair_compact_gpu.cu \
	pair_sort_gpu.cu \
	bonded_migration.fpp \
	bonded_migration_gpu.cu \
	bonded_residency_gpu.cu \
	exclusion_gpu.fpp \
	exclusion_gpu_ops.cu \
	pair_prepare.fpp \
	pair_prepare_gpu.cu \
	pme_mesh_sum.fpp \
	pme_host_pack.fpp \
	gpu_common.h \
	spdyn.fpp
else
spdyn_SOURCES += \
	spdyn.fpp
endif

spdyn_ADDLIB = ../lib/lib.a

SUFFIXES = .fpp

EXTRA_DIST = \
	Makefile.depends 

DEPENDS	= Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(spdyn_OBJECTS) *.f90 *~ *.mod

sp_enefunc_charmm.o: domain_masks.o
exclusion_oracle.o: sp_enefunc_str.o
pair_sparse.o: exclusion_reset.o
sp_enefunc_charmm.o: exclusion_reset.o exclusion_oracle.o

.fpp.o:
if USE_NVHPC
	$(FPP2F90_NVHPC) $*.fpp > $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP2F90) $*.f90
	mv $*.cpp.f90 $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
	$(FPP2F90) $*.fpp $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
endif
endif

if USEKCOMP
sp_pairlist_fugaku.o: FCFLAGS=$(FCFLAGS_KERNELS_PAIRLIST)
sp_energy_dihedrals.o: FCFLAGS=$(FCFLAGS_DIHEDRAL)
sp_energy_nonbond_fugaku.o: FCFLAGS=$(FCFLAGS_KERNELS)
endif
if !USEKCOMP
sp_energy_pme_noopt_1dalltoall.o sp_energy_pme_noopt_2dalltoall.o sp_energy_pme_noopt_slab.o: FCFLAGS=$(FCFLAGS_NOOPT)
endif

if CUDAGPU
.cu.o:
	$(NVCC) $(NVCCFLAG) --ptxas-options=-v -lineinfo $*.cu
endif

sp_pairlist_gpu.o: pair_host.o pair_sparse.o
sp_enefunc_charmm.o: pair_sparse.o

spdyn$(EXEEXT) : $(spdyn_ADDLIB) $(spdyn_OBJECTS)
	$(FC) -o spdyn $(spdyn_OBJECTS) \
		$(spdyn_ADDLIB) $(LDFLAGS) $(PME_CUFFT_LIBS) $(STEP_CUDA_LIBS)

$(spdyn_ADDLIB):
	cd ../lib; ${MAKE}

depend: clean_depend
	python3 ../../fortdep.py -e ../lib -f ".fpp" *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

sp_update_domain.o sp_enefunc.o sp_md_respa.o: maintenance_probe.o
sp_energy.o sp_energy_nonbond_gpu.o: maintenance_probe.o
if CUDAGPU
gpu_step.o: maintenance_probe.o
endif

if CUDAGPU
translation_gpu.o gpu_force_assembly.o water_migration.o hgr_migration.o: precision_types.h
gpu_pme.o long_force_gpu.o: precision_types.h
gpu_state_ops.o observables_gpu.o maintenance_state_gpu.o metadata_import_gpu.o: precision_types.h
water_migration.o hgr_migration.o atom_outgoing_gpu.o incoming_gpu_ops.o: migration_coordinates.h
force_topology_gpu.o bonded_migration_gpu.o bonded_residency_gpu.o bonded_incoming_hbm.o: precision_types.h
gpu_sp_energy.o nonbond_mirror_gpu.o nonbond_virial_gpu.o: precision_types.h
gpu_coord_pack.o remap_coordinate_gpu.o state_cycle_gpu.o: precision_types.h
incoming_gpu_ops.o atom_outgoing_gpu.o migration_resident_gpu.o pair_prepare_gpu.o exclusion_gpu_ops.o: precision_types.h
incoming_gpu_ops.o atom_outgoing_gpu.o migration_resident_gpu.o pair_prepare_gpu.o exclusion_gpu_ops.o improper_outgoing_gpu.o: migration_views.h
bonded_incoming_gpu.o bonded_incoming_hbm.o bonded_incoming_cooperative.o improper_outgoing_gpu.o force_topology_gpu.o: bonded_incoming_view.h
force_topology_bindings.o: FCFLAGS += -cuda
bonded_force.o: force_topology_bindings.o
metadata_import.o: FCFLAGS += -cuda
metadata_import.o: sp_domain_str.o sp_constraints_str.o
gpu_step.o: metadata_import.o
metadata_import_gpu.o: migration_views.h
state_cycle.o: FCFLAGS += -cuda
state_cycle.o: sp_domain_str.o
gpu_step.o: state_cycle.o
state_cycle_gpu.o: migration_views.h
observables.o: FCFLAGS += -cuda
observables.o: sp_domain_str.o
gpu_step.o gpu_thermo.o sp_dynvars.o: observables.o
improper_outgoing.o: sp_domain_str.o sp_enefunc_str.o bonded_incoming.o
sp_migration.o: improper_outgoing.o
pair_prepare.o: sp_domain_str.o
sp_pairlist_gpu.o: pair_prepare.o

exclusion_gpu.o: sp_domain_str.o sp_enefunc_str.o sp_constraints_str.o domain_masks.o
sp_enefunc_charmm.o: exclusion_gpu.o
endif

incoming_oracle.o: sp_domain_str.o sp_constraints_str.o
sp_migration.o: incoming_oracle.o

sp_pairlist_gpu.o: pair_sort.o

if CUDAGPU
bonded_migration.o: sp_domain_str.o sp_enefunc_str.o
sp_migration.o: bonded_migration.o
endif

-include $(DEPENDS)

if CUDAGPU
gpu_sp_energy.o nonbond_inter_staged.o: nonbond_inter_staged.h precision_types.h gpu_common.h
endif
