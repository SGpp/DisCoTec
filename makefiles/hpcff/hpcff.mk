###########################################################################
### architecture dependent GENE makefiles, HPC-FF version               ###
### http://www2.fz-juelich.de/jsc/juropa/                               ###
###########################################################################

###########################################################################
### make sure you have the following line in your ~/.bashrc script      ###
### module unload intel mkl                                             ###
### module load intel/12.1.4 mkl                                        ###
### module load PETSc/3.4.2-basic-complex_opt_intel12.1.4               ###
###########################################################################

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
#has to be switched on
SCALAPACK = yes
DEBUG = no
USE_PERFLIB = none
FUTILS = yes
PRODRUN = no

MKLVERSION = 10.2

#memory per core
MB_PER_CORE=2800
FCK_PREPROC = 'MP_PER_CORE=2800'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

CHIP = Nehalem

MKLROOT = $(MKLINCLUDE)/../

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

ifeq ($(FFTLIB),fftw)
 ifeq ($(PRECISION),double)
   FFTWPATH = /usr/local/fftw/3.2.2
   LIBS += -L$(FFTWPATH)/lib -lfftw3
 else
   FFTWPATH = /usr/local/fftw/3.2.2_float
   LIBS += -L$(FFTWPATH)/lib -lfftw3f
 endif
 INCPATHS += -I$(FFTWPATH)/include
endif


###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 ifeq ($(PRECISION),double)
#  PETSC_DIR = /lustre/jhome6/fsngeul/fsgeul00/soft/petsc
  PETSC_ARCH = arch-installed-petsc
  SLEPC_DIR = /lustre/jhome6/fsngeul/fsgeul00/soft/slepc-3.4.3/
  include $(SLEPC_DIR)/conf/slepc_common
  INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
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
 HDF5PATH = /usr/local/hdf5/v1.8.5/
 HDF5_LIBS =  -L$(HDF5PATH)/lib  -lhdf5_fortran -lhdf5 -lz /usr/local/szip/lib/libsz.a
 LIBS += -L$(FUTILSDIR) -lfutils $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)

endif

ifeq ($(USE_PERFLIB),perf)
 PREPROC += -DWITHPERF
 LIBS += -lperf_r
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


###########################################################################
### Running                                                             ###
###########################################################################

MPRUN = mpiexec -env MKL_SERIAL yes -env OMP_NUM_THREADS $(OMP_NUM_THREADS) -np=$(N_PES) ./$(EXEC)

###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
#NOOPTLIST = init_cond.o