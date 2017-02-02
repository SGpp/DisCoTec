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

MPI_IO = no
SLEPC = no
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = yes
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = no
MKLVERSION = 10.3
USE_MEMORY_CHECKER=none
USE_C_NONLIN=no
USE_CUDA_NONLIN=no

#definition of the platform
# memory per core
MB_PER_CORE=5000
# chip type
CHIP = Nehalem

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
MKLROOT ? =$(MKL_HOME)

#include $(BASEDIR)/makefiles/compilers/$(COMPILER).def


### FFT LIBRARY
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

ifeq ($(USE_PERFLIB),perf)
 LIBS += $(PERFLIB_TS)
 PREPROC += -DWITHPERF=1
endif

ifeq ($(DEBUG),no)
# this setting is important, otherwise is the host code of the
# .cu files unoptimized by gcc
  CUDAFLAGS += -O3
endif

ifeq ($(INCLUDE_SYMBOLS),yes)
	CUDAFLAGS += -g -G
endif

ifeq ($(OPENMP),yes)
	PREPROC += -DWITHOMP_BLOCKLOOP
endif

###########################################################################
### Linking                                                             ###
###########################################################################
MPLD = $(MPFC)
LD = $(FC)

NOOPTLIST= 

NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

#forcheck:
#	forchk  -I $(FCK_INCPATHS) -l mylistfile \
#		-ff -i4 -dp -allc -decl -nshsrc -ancmpl -anref -shmod \
#		$(F90PRENAMES) $(PPDIR)/gene.f90 -lib $(FCKDIR)/share/forcheck/MPI.flb
