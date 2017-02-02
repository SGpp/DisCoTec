###########################################################################
### architecture dependent GENE makefiles, gp cluster version          ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load impi/4.0.3 intel/12.1 mkl/10.3 fftw petsc slepc				  # 
# module load hdf5-mpi perflib/2.0                                               #
#								          #
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
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = yes
MKLVERSION = 10.3
USE_MEMORY_CHECKER=none
USE_C_NONLIN=no
USE_CUDA_NONLIN=no


#memory per core
PREPROC = -D'MB_PER_CORE=7000'
FCK_PREPROC = 'MP_PER_CORE=1900'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
CHIP = SandyBridge
MKLROOT ?= $(MKL_HOME)
OPTLEVEL=3

#include $(BASEDIR)/makefiles/compilers/$(COMPILER).def

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   FCK_INCPATHS := $(FCK_INCPATHS):$(FFTW_HOME)/include
 ifeq ($(PRECISION),double)
   LIBS += $(FFTW_HOME)/lib/libfftw3.a
 else
   LIBS += $(FFTW_HOME)/lib/libfftw3f.a
 endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
  FCK_PREPROC := $(FCK_PREPROC):WITHSCAL
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
 INCPATHS += -I$(FUTILSDIR)
endif

MEMCHECK_OBJ=
ifeq ($(USE_MEMORY_CHECKER),own)
	MEMCHECK_OBJ = $(OBJDIR)/set_malloc_hooks.o
	PREPROC += -DMEMORY_CHECKER
endif

###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################
#PRECOMILER SWITCHES

ifeq ($(USE_C_NONLIN),yes)
	PREPROC += -DWITH_C_NONLIN
endif
ifeq ($(USE_CUDA_NONLIN),yes)
	PREPROC += -DWITH_CUDA_NONLIN
	LIBS += -L$(CUDA_HOME)/lib64 -lcudart -lcufft
endif

NVCC = nvcc
CUDAFLAGS += -gencode arch=compute_20,code=sm_20
CUDAFLAGS += -g -G

ifeq ($(USE_PERFLIB),perf)
 LIBS += -L$(PERFLIB_HOME)/lib -looperf -Wl,-rpath,$(PERFLIB_HOME)/lib
 PREPROC += -DWITHPERF=1
 INCPATHS += -I$(PERFLIB_HOME)/include
endif

ARCHIVE = ar r

ifeq ($(DEBUG),no)
 CUDAFLAGS += -O3
endif

ifeq ($(COMPILER_REPORTS),yes)
#	FFLAGS += $(REPORT_FFLAGS)
endif

ifeq ($(OPENMP),yes)
	FFLAGS += $(OPENMP_FFLAGS)
	PREPROC += -DWITHOMP -DWITHOMP_NONLIN
endif


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

ifeq ($(FFTLIB),mkl)
 $(OBJDIR)/mkl_dfti.o: $(MKLROOT)/include/mkl_dfti.f90
endif

NOOPTLIST= 

NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

