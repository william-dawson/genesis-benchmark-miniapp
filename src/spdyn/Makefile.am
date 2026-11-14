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
gpu_thermo.o: sp_constraints.o sp_domain_str.o sp_constraints_str.o
water_bindings.o: FCFLAGS += -cuda
water_bindings.o: sp_domain_str.o sp_boundary_str.o
gpu_step.o sp_migration.o: water_bindings.o
hgr_bindings.o: FCFLAGS += -cuda
hgr_bindings.o: sp_domain_str.o sp_boundary_str.o sp_constraints_str.o
gpu_step.o sp_migration.o: hgr_bindings.o
coord_bindings.o: FCFLAGS += -cuda
coord_bindings.o: sp_domain_str.o
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
	gpu_state_ops.cu \
	gpu_force_assembly.cu \
	coord_bindings.fpp \
	gpu_coord_pack.cu \
	hgr_bindings.fpp \
	hgr_migration.cu \
	water_bindings.fpp \
	water_migration.cu \
	gpu_thermo.fpp \
	gpu_step.fpp \
	gpu_sp_energy.cu \
	gpu_pme.cu \
	pair_sparse_gpu.cu \
	pair_compact_gpu.cu \
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

incoming_oracle.o: sp_domain_str.o sp_constraints_str.o
sp_migration.o: incoming_oracle.o

-include $(DEPENDS)
