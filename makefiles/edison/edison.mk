###########################################################################
### architecture dependent GENE makefile for EDISON (NERSC)             ###
###  http://www.nersc.gov/systems/edison-cray-xc30/                     ###
###########################################################################

###########################################################################
### List of required environment settings/modules:                      ###
###                                                                     ###
### module load fftw cray-libsci cray-petsc-complex                     ###
### module load gv gnuplot idl						###
### module load hdf5-parallel                                           ###
### (copy lines, e.g., in your ~/.cshrc.ext)                            ###
###                                                                     ###
### slepc is currently installed and accessible on the home folder of   ###
### one of the Edison GENE users as a temporary workaround. NERSC       ###
### installation of the slepc module on Edison has been requested       ###
###									###	
### The hdf5-parallel module is marked as deprecated; the alternative   ###
### cray-hdf5-parallel module, however, doesn't seem to work with intel ###
### compilers at the moment                                             ###  
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use cray wrappers for compiler names:
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP = Xeon_E5_IvyBridge

ARCHIVE = ar r

NTASKS ?= 24
MPRUN = aprun -n $(N_PES) -N $(NTASKS) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = fftw

PRECISION = double

DEBUG= no

SLEPC= yes

SCALAPACK = yes

OPENMP = no

USE_PERFLIB = none

DIAG_MPI_IO= no

FUTILS = no

MB_PER_CORE=2400
#officially: 2.67GB/core
#however, the MPI environment itself typically consumes some fraction 

INCLUDE_SYMBOLS = no

COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = 

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = 

#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE = 
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 ifeq ($(FFTW_INC),)
  $(error run module load fftw first)
 endif
 INCPATHS += -I$(FFTW_INC)
 LIBS += -L$(FFTW_DIR) -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# set by module load petsc-complex:
 PETSC_DIR ?= 
 SLEPC_DIR ?= /global/u1/t/tbg/soft/edison/slepc-3.4.3
 SLEPC_DIR ?= /global/u1/j/jcitrin/gene11_release16/slepc-3.3-p3
 ifeq ($(strip $(PETSC_ARCH)),)
  PETSC_ARCH = arch-installed-petsc
 endif

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS += 
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = 
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST = 
