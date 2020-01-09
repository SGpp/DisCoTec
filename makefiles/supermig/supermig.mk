###########################################################################
### architecture dependent GENE makefiles, BOB cluster version          ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load prace							  #
# module load mpi.intel fftw petsc/3.1c slepc				  # 
# module load hdf5-mpi perf                                               #
# 
# IBM MPI:
# module load poe/1.2 mpi.ibm/1.2
#								          #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=intel

# comment the following two lines for use of Intel MPI
MPFC = mpif90
MPCC = mpicc

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw
PRECISION= double
OPENMP = no

MPI_IO = no
SLEPC = yes
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = yes
MKLVERSION = 10.3

#memory per core
MB_PER_CORE=1900
FCK_PREPROC = 'MP_PER_CORE=1900'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
CHIP = Nehalem
MKLROOT = $(MKL_BASE)

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += $(FFTW_INC)
   LIBS += $(FFTW_LIB)
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
#   INCPATHS += $(SLEPC_INCLUDE) $(PETSC_INC) -I$(PETSC_DIR)/$(PETSC_ARCH)/include
   LIBS += $(SLEPC_LIB)
   PREPROC += -DWITHSLEPC
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
  PREPROC += -DWITHSCAL
  FCK_PREPROC := $(FCK_PREPROC):WITHSCAL
else
#Insert only BLAS library
  LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)
endif


###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################

ifeq ($(USE_PERFLIB),perf)
 LIBS += $(PERFLIB_TS)
 PREPROC += -DWITHPERF=1
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

#LDFLAGS += -Xlinker -M
#LIBS += -L/lrz/sys/intel/ifort_121_339/lib/intel64 -lsvml

###########################################################################
### Running                                                             ###
###########################################################################
MPRUN = export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=yes;\
	mpiexec -n $(N_PES) ./$(EXEC)

###########################################################################


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST= 

preproc: $(F90PRENAMES) $(PPDIR)/gene.f90
