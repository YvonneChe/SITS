#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999, 2004, 2008 **
#                Regents of the University of California               **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************

# Makefile for Amber 14  PMEMD 
SHELL=/bin/sh
AMBER=Amber14
CONFIG_FILE=../../config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA
CUDA_DPFP_TOKEN=Duse_DPFP
CUDA_SPFP_TOKEN=Duse_SPFP
CUDA_SPXP_TOKEN=Duse_SPXP
MIC_TOKEN=mmic
MIC_OFFLOAD_TOKEN=DMIC_offload

# Platform-specific info should be found in ../../config.h
include $(CONFIG_FILE)

OBJS=   gbl_constants$(OSFX) gbl_datatypes$(OSFX) state_info$(OSFX) file_io_dat$(OSFX) \
        mdin_ctrl_dat$(OSFX) mdin_emil_dat$(OSFX) mdin_ewald_dat$(OSFX) mdin_subunit_dat$(OSFX) mdin_its_dat$(OSFX) mdin_debugf_dat$(OSFX) prmtop_dat$(OSFX) \
        inpcrd_dat$(OSFX) dynamics_dat$(OSFX) emil$(OSFX) img$(OSFX) nbips$(OSFX) offload_allocation$(OSFX) \
        parallel_dat$(OSFX) parallel$(OSFX) gb_parallel$(OSFX) \
        pme_direct$(OSFX) pme_recip_dat$(OSFX) pme_slab_recip$(OSFX) pme_blk_recip$(OSFX) \
        pme_slab_fft$(OSFX) pme_blk_fft$(OSFX) pme_fft_dat$(OSFX) fft1d$(OSFX) \
        bspline$(OSFX) pme_force$(OSFX) pbc$(OSFX) nb_pairlist$(OSFX) \
        nb_exclusions$(OSFX) cit$(OSFX) dynamics$(OSFX) bonds$(OSFX) angles$(OSFX) dihedrals$(OSFX) \
        extra_pnts_nb14$(OSFX) runmd$(OSFX) loadbal$(OSFX) shake$(OSFX) prfs$(OSFX) mol_list$(OSFX) \
        runmin$(OSFX) constraints$(OSFX) axis_optimize$(OSFX) gb_ene$(OSFX) veclib$(OSFX) gb_force$(OSFX) its_force$(OSFX) \
        timers$(OSFX) pmemd_lib$(OSFX) runfiles$(OSFX) file_io$(OSFX) \
        AmberNetcdf$(OSFX) bintraj$(OSFX) binrestart$(OSFX) pmemd_clib$(OSFX) \
        pmemd$(OSFX) random$(OSFX) degcnt$(OSFX) erfcfun$(OSFX) nmr_calls$(OSFX) nmr_lib$(OSFX) \
        get_cmdline$(OSFX) master_setup$(OSFX) pme_alltasks_setup$(OSFX) pme_setup$(OSFX) \
        ene_frc_splines$(OSFX) gb_alltasks_setup$(OSFX) nextprmtop_section$(OSFX) \
        angles_ub$(OSFX) dihedrals_imp$(OSFX) cmap$(OSFX) charmm$(OSFX) charmm_gold$(OSFX) \
        findmask$(OSFX) remd$(OSFX) multipmemd$(OSFX) remd_exchg$(OSFX) amd$(OSFX) ti$(OSFX) gbsa$(OSFX) \
        barostats$(OSFX) scaledMD$(OSFX) constantph$(OSFX) energy_records$(OSFX) constantph_dat$(OSFX) \
        relaxmd$(OSFX) sgld$(OSFX) emap$(OSFX)  

install: configured_serial $(BINDIR)/pmemd.SITS$(SFX)

parallel: configured_parallel $(BINDIR)/pmemd.SITS.MPI$(SFX)
	@echo " "

cuda: configured_cuda pmemd.SITS.cuda$(SFX)
	@( \
	  if grep $(CUDA_SPXP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.SITS.cuda$(SFX) $(BINDIR)/pmemd.SITS.cuda_SPXP$(SFX) ;\
	  elif grep $(CUDA_DPFP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.SITS.cuda$(SFX) $(BINDIR)/pmemd.SITS.cuda_DPFP$(SFX) ;\
	  else \
	    mv pmemd.SITS.cuda$(SFX) $(BINDIR)/pmemd.SITS.cuda_SPFP$(SFX) ;\
	    cd $(BINDIR) ; ln -f -s pmemd.SITS.cuda_SPFP$(SFX) pmemd.SITS.cuda$(SFX);\
	  fi ;\
	)

cuda_parallel: configured_cuda configured_parallel pmemd.SITS.cuda.MPI$(SFX)
	@( \
	  if grep $(CUDA_SPXP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.SITS.cuda.MPI$(SFX) $(BINDIR)/pmemd.SITS.cuda_SPXP.MPI$(SFX) ;\
	  elif grep $(CUDA_DPFP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.SITS.cuda.MPI$(SFX) $(BINDIR)/pmemd.SITS.cuda_DPFP.MPI$(SFX) ;\
	  else \
	    mv pmemd.SITS.cuda.MPI$(SFX) $(BINDIR)/pmemd.SITS.cuda_SPFP.MPI$(SFX) ;\
	    cd $(BINDIR) ; ln -f -s pmemd.SITS.cuda_SPFP.MPI$(SFX) pmemd.SITS.cuda.MPI$(SFX);\
	  fi ;\
	)

mic: configured_mic $(BINDIR)/pmemd.mic_native$(SFX)

mic_parallel: configured_mic configured_parallel $(BINDIR)/pmemd.mic_native.MPI$(SFX)

mic_offload: configured_mic_offload $(BINDIR)/pmemd.mic_offload.MPI$(SFX)

$(BINDIR)/pmemd.SITS$(SFX): $(OBJS) $(EMIL) 
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

$(BINDIR)/pmemd.SITS.MPI$(SFX): $(OBJS) $(EMIL) 
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

$(BINDIR)/pmemd.mic_native$(SFX): $(OBJS) $(EMIL)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) 

$(BINDIR)/pmemd.mic_offload.MPI$(SFX): $(OBJS) $(EMIL)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

$(BINDIR)/pmemd.mic_native.MPI$(SFX): $(OBJS) $(EMIL)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

pmemd.SITS.cuda$(SFX): $(OBJS) $(PMEMD_CU_LIBS) $(EMIL)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
      $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

pmemd.SITS.cuda.MPI$(SFX): $(OBJS) $(PMEMD_CU_LIBS) $(EMIL)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
     $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF)

EMIL:
	$(MAKE) -C ../../../AmberTools/src/emil install

depends:
	../makef90depends > f90depends

clean:
	$(RM) -f *.o *.obj *.mod pmemd.SITS.MPI$(SFX) pmemd.SITS$(SFX) pmemd.SITS.cuda$(SFX) pmemd.SITS.cuda.MPI$(SFX) pmemd.mic_native$(SFX) pmemd.mic_offload.MPI$(SFX) pmemd.mic_native.MPI$(SFX) *.d work.pc* 
	$(MAKE) -C ./cuda clean

# Control the suffixes available; this essentially eliminates built-in
# inference rules, done here to improve portability.
.SUFFIXES:
.SUFFIXES: .F90 .c $(OSFX)

.F90$(OSFX):
	$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.F90

.c$(OSFX):
	$(PMEMD_CC) $(PMEMD_COPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.c

$(PMEMD_CU_LIBS):
	$(MAKE) -C ./cuda

AmberNetcdf$(OSFX): ../../../AmberTools/src/lib/AmberNetcdf.F90
	$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c -o AmberNetcdf$(OSFX) ../../../AmberTools/src/lib/AmberNetcdf.F90

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE)!" ;\
	    exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type parallel!" ;\
		echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
		exit 2 ;\
	} ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -mpi." ;\
		exit 2 ;\
	fi ;\
	)
	@(if grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type cuda, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -cuda." ;\
		exit 2 ;\
	fi ;\
	)
	@(if grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type mic_native, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -mic_native." ;\
		exit 2 ;\
	fi ;\
	)
	@(if grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type mic_offload, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -mic_offload." ;\
		exit 2 ;\
	fi ;\
	)

configured_cuda: configured
	@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type cuda!" ;\
		echo "  Rerun $(CONFIG_COMMAND) and specify a cuda implementation." ;\
		exit 2 ;\
		} ;\
	)

configured_mic: configured
	@(grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic_native!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic_native flag." ;\
		exit 2 ;\
		} ;\
	)

configured_mic_offload: configured
	@(grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic_offload!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic_offload flag." ;\
		exit 2 ;\
		} ;\
	)

## Dependencies for f90 modules and include files generated by makef90depends:

include f90depends

