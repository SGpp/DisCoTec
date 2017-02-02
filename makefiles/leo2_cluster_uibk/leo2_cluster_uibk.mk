###########################################################################
### architecture dependent GENE makefile for the LEO-II cluster at UIBK ###
### see http://www.uibk.ac.at/zid/systeme/hpc-systeme/leo2/ for details ###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= no 

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = no

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

#performance profiling library 
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required for some geometry interfaces
FUTILS = no

#memory per core
PREPROC= -D'MB_PER_CORE=750'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR) -I$(UIBK_INTEL_MKL_INC)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = -L$(UIBK_INTEL_MKL_LIB) -lmkl_lapack \
	-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

#FFT LIBRARY
#specify at least one of the following choices
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
 INCPATHS += -I$(UIBK_FFTW_INC)
 LIBS +=  -L$(UIBK_FFTW_LIB) -lfftw3 
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if those variables are not
# set by modules (make sure you use the complex versions)
# PETSC_ARCH = 
# PETSC_DIR = 
# SLEPC_DIR = 

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS +=
  PREPROC += -DWITHSCAL
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = 
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)

endif


###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpif90

#ARCHIVE command
#Note: On IBM machines you might need to add -X$(BITMODE) to the "ar" call
ARCHIVE = ar r

#FORTAN COMPILER FLAGS
FFLAGS = 

#PRECOMILER SWITCHES
PREPROC += -DMPI $(SVNDEF)

#Uncomment the following line if your FORTRAN compiler does not provide
#erf(x) as intrinsic function:
#PREPROC += -DEXTERNAL_ERF

#LINKER (usually same as MPFC)
MPLD =$(MPFC)

#LDFLAGS
LDFLAGS =

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
ifeq ($(DEBUG),yes)
#switch off any optimization 
 OPT=	-O0
#***TO DO***: set some flags for backtracing, e.g.:
#OPT+= -g -C etc
else
 OPT= 
#***TO DO***: check the compiler documentation for optimization flags
endif

FFLAGS += $(OPT)
LDFLAGS += $(OPT)

#OPENMP
ifeq ($(OPENMP),yes)
FFLAGS +=
LDFLAGS += 
PREPROC +=-DWITHOMP
endif

ifeq ($(PRECISION),double)
FFLAGS += -r8
PREPROC +=-DDOUBLE_PREC
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif

###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
	mpirun -np $(N_PES) ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#additional compilations, e.g., of MKL library files, might be 
#necessary on some machines (see for instance bob_cluster)
#
#Furthermore, optimization options might be adapted for certain files here
