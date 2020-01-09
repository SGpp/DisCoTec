###########################################################################
### Machine dependent GENE makefile, Lindgren Cray XT6m version         ###
### (http://www.pdc.kth.se/resources/computers/lindgren)                ###
###########################################################################
#                                                                         #
# export MACHINE=lindgren                                                 #
#                                                                         #
## Load the following modules (either by hand or in your .bashrc/.tcshrc) #
#                                                                         #
# module unload PrgEnv-pgi petsc-complex slepc-complex                    #
# module load PrgEnv-intel petsc-complex/3.3.00                           #
#                                                                         #
## and maybe								  #
# module load xt-libsci                                                   #
#                                                                         #
# slepc is currently installed and accessible on the home folder of       #
# one of the lindgren GENE users as a temporary workaround.               #
#                                                                         #  
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#use cray wrappers for compiler names:
MPFC = ftn
MPCC = cc
COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)
#CHIP = AMD_Opteron_6100

ARCHIVE = ar r

#NTASKS ?= 24
MPRUN = OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
        pwd;\
        aprun -n $(N_PES) ./$(EXEC)

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

MB_PER_CORE=1100
#officially: 1.33GB/core
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
 SLEPC_DIR = /afs/pdc.kth.se/home/t/tegnered/Public/slepc-3.3-p4/

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
