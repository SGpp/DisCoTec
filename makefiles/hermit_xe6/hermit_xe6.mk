###########################################################################
### architecture dependent GENE makefile for HERMIT                     ###
### https://wickie.hlrs.de/platforms/index.php/Cray_XE6                 ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
### switch to desired PrgEnv,e.g.:                                      ###
### module unload PrgEnv-cray && module load PrgEnv-intel               ###
### then                                                                ###
### module load fftw                                                    ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use Hermit wrappers
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

CHIP = interlagos
#note: optimization flags are also set by the compiler wrappers themselves;
#hence, there's no urgent need (in theory) to apply further corresponding 
#flags in our compiler makefile

ARCHIVE = ar r

NTASKS ?= 32
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

MB_PER_CORE=900

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
 PETSC_ARCH ?= $(COMPILER)-double-complex
 PETSC_DIR ?= /zhome/academic/HLRS/pri/iprtgoer/soft/petsc
 SLEPC_DIR ?= /zhome/academic/HLRS/pri/iprtgoer/soft/slepc

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS += -lsci_$(COMPILER)
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
