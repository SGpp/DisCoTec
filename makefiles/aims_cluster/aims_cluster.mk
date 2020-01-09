###########################################################################
### architecture dependent GENE makefiles, AIMS cluster version         ###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

#COMPILER: INTEL or LAHEY
COMPILER = intel

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw
PRECISION= double
OPENMP = no

SLEPC = no
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = yes
MKLVERSION = 10.3

################## some libraries are not yet available for  Lahey #######
ifeq ($(COMPTYPE),LAHEY)
 FFTLIB = fftw
 SLEPC  = no
 FUTILS  = no
endif

#memory per core
PREPROC= -D'MB_PER_CORE=2200'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################
CHIP = 
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
#Insert BLAS library
#LIBS += 

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
  ifeq ($(COMPTYPE),LAHEY)
    FFTWPATH = /afs/ipp/home/t/tisd/bob_inst/lahey
    INCPATHS += -I$(FFTWPATH)/include
    LIBS += -L$(FFTWPATH)/lib
  endif
  ifeq ($(COMPTYPE),INTEL)
    FFTWPATH = $(FFTW_HOME)
    INCPATHS += -I$(FFTWPATH)/include
    LIBS += -L$(FFTWPATH)/lib
  endif
  ifeq ($(PRECISION),double)
   LIBS += -lfftw3
  else
   LIBS += -lfftw3f
  endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
#   PETSC_ARCH = amd64_sles11-mpicc_11.1-mpif90_11.1-64-mkl-double-complex.petsc-3.1-p4
#   PETSC_DIR = ${PETSC_DIR}
#   SLEPC_DIR = ${SLEPC_DIR}
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
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
 HDF5PATH = /afs/ipp/home/t/tbg/soft/@sys/hdf5
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)

endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
ifeq ($(USE_PERFLIB),perf)
 LIBS += -lperf_r
endif
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
MPRUN = export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=yes;\
	mpiexec -n $(N_PES) ./$(EXEC)

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

NOOPTLIST= 

preproc: $(F90PRENAMES) $(PPDIR)/gene.f90
