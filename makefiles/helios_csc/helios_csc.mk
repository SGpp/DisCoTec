###########################################################################
### architecture dependent GENE makefiles for BULL supercomputer helios ###
### https://www.iferc-csc.org/                                          ###
###########################################################################
#                                                                         #
# the following modules have to be loaded first (e.g., in ~/.bash_profile)#
# module load bullxmpi intel                                              # 
# module load fftw hdf5_p/1.8.8                                           #
## module load petsc/3.2-p6/complex-avx slepc/3.2-p3/complex-avx          #
# module load mxml adios                                                  #
#                                                                         #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER = intel
CHIP = SandyBridge

#ARCHIVE command
ARCHIVE = ar r

MPRUN = export OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
        export MKL_SERIAL=yes;\
        mpiexec -n $(N_PES) $(EXEC)


# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

PRECISION= double

# Choose MPIVENDOR=intel for Intel MPI and bullx for Bull MPI
ifneq ($(I_MPI_ROOT),)
MPIVENDOR = intelmpi-4.0.3
else
MPIVENDOR = bullxmpi-1.1.11.1
endif

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= yes

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

#performance profiling library 
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required for some geometry interfaces
FUTILS = no

# HLST ADIOS CHECKPOINT
ifneq ($(ADIOS_DIR),)
 HAC = yes
else
 HAC = no
endif

#Provide an upper (RAM) memory limit for the code internal performance 
#optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=3500
#Note: the system itself often requires a significant fraction 
#      which should be taken into account

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

MKLPATH = /csc/softs/intel/mkl/
MKLROOT = $(MKLPATH)
MKLLIBPATH = $(MKLPATH)/lib/intel64

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = 

#FFT LIBRARY

ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(FFTW_DIR)/include
 LIBS += -L$(FFTW_DIR)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH = arch-linux2-cxx-opt
 PETSC_DIR = /csc/home3/tbg/soft/petsc/
 SLEPC_DIR = /csc/home3/tbg/soft/slepc/

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

ifeq ($(SCALAPACK),yes)
 LIBS += $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
                -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
                $(MKLROOT)/lib/intel64/libmkl_sequential.a \
                $(MKLROOT)/lib/intel64/libmkl_core.a
 ifeq ($(findstring intel,$(MPIVENDOR)),intel)
	LIBS += $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a
 else
	LIBS += $(MKLROOT)/lib/intel64/libmkl_blacs_openmpi_lp64.a 
 endif
 LIBS += -Wl,--end-group -lpthread -lm
else
#Insert BLAS library
LIBS += -L$(MKLLIBPATH)
LIBS += -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a \
$(MKLLIBPATH)/libmkl_sequential.a $(MKLLIBPATH)/libmkl_core.a \
-Wl,--end-group -lpthread
endif

ifeq ($(findstring intel,$(MPIVENDOR)),intel)
#hdf5 is not installed yet for intelmpi
FUTILS=no
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 
# HDF5_LIBS += -lhdf5hl_fortran -lhdf5_hl -lz

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR) -I$(HDF5PATH)/include

endif

ifeq ($(HAC),yes)
 #  HLST-ADIOS-CHECKPOINT interface
 ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
 ADIOS_LIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
endif


ifeq ($(USE_PERFLIB),perf)
 LIBS += -L/csc/home0/dannert/helios_inst/lib -looperf_r -lpfm -lstdc++ 
 LD_LIBRARY_PATH += /csc/home0/dannert/helios_inst/lib
endif


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################


