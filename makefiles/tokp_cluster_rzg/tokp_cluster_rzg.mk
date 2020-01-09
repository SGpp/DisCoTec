###########################################################################
### architecture dependent GENE makefiles, TOK-P cluster version        ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load intel mkl fftw petsc slepc                                  #
# module load hdf5-mpi                                                    #
# 									  #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=intel

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw
PRECISION= double
OPENMP = no

SLEPC = yes
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = yes
PRODRUN = yes
MKLVERSION = 10.3

#memory per core
MB_PER_CORE=3500
FCK_PREPROC = 'MP_PER_CORE=3500'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
CHIP = SandyBridge
ifneq ($(MKL_BASE),)
	MKLROOT ?= $(MKL_BASE)
else
ifneq ($(MKL_HOME),)
	MKLROOT ?= $(MKL_HOME)
endif
endif

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS += 
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   LIBS += -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   PETSC_ARCH = 
   PETSC_DIR = /tokp/work/tbg/soft/petsc/
   SLEPC_DIR = /tokp/work/tbg/soft/slepc/
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
#  compiling with SLEPc requires libsvml to avoid unresolved dependencies
   LIBS += -L$(IFORT_BASE)/compiler/lib/intel64 -lsvml
   LIBS += $(SLEPC_LIB)
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
else
#Insert only BLAS library
  LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR) -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)
endif


###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################

ifeq ($(USE_PERFLIB),perf)
 LIBS += -L$(PERFLIB_HOME)/lib -looperf -lpfm
 PREPROC += -DWITHPERF=1
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

#LDFLAGS += -Xlinker -M

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
