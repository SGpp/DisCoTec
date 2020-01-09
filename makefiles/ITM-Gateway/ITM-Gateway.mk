###########################################################################
### architecture dependent GENE makefile - ITM Gateway version          ###
### https://itm.ipp.mpg.de/						###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
###                                                                     ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

COMPILER = intel
CHIP = SandyBridge

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = mpiexec -n $(N_PES) ./$(EXEC)
#MPRUN = mpirun -np $(N_PES) ./$(EXEC)


###########################################################################
### SWITCHES                                                            ###
###########################################################################


FFTLIB = fftw
PRECISION= double
DEBUG= no
SLEPC= yes
SCALAPACK = no
OPENMP = no
USE_PERFLIB = none
FUTILS = yes
HAC = yes

MB_PER_CORE=3500
#officially, 4GB/core


INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

OBJECTCODE=amd64_intel_12
ITM_INCLUDE += $(shell eval-pkg-config --cflags itmtypes-${OBJECTCODE})
ITM_LIBS = $(shell eval-pkg-config --libs itmtypes-${OBJECTCODE})
ITM_INCLUDE += $(shell eval-pkg-config --cflags itmconstants-${OBJECTCODE})
ITM_LIBS += $(shell eval-pkg-config --libs itmconstants-${OBJECTCODE})
ITM_INCLUDE += $(shell eval-pkg-config --cflags xmllib-${OBJECTCODE})
ITM_LIBS += $(shell eval-pkg-config --libs xmllib-${OBJECTCODE})
ITM_INCLUDE += $(shell eval-pkg-config --cflags ual-${OBJECTCODE})
ITM_LIBS += $(shell eval-pkg-config --libs ual-${OBJECTCODE})


#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = $(ITM_INCLUDE)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = -L/afs/ipp/itm/switm/lib/lapack/3.4.2/lib/amd64_intel_12 -llapack -lrefblas
LIBS+=$(ITM_LIBS)


#FFT LIBRARY
#fill at least one of the following choices: mkl, essl, fftw
#and specify this choice as FFTLIB in the header

ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE = 
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),fftw)
 FFTW_HOME =/afs/rzg/common/soft/fftw/fftw-3.3.3/amd64_sles11/intel-12.1/impi-4.1/
 INCPATHS += -I$(FFTW_HOME)/include
ifeq ($(PRECISION),double)
   LIBS += -L$(FFTW_HOME)/lib/ -lfftw3
 else
   LIBS += -L$(FFTW_HOME)/lib/ -lfftw3f
 endif
 LIBS += -lm
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# uncomment and fill the following if those variables are not
# set by modules (make sure, you use the complex versions)
PETSC_ARCH ?= intel12.1-impi4.1.0-mkl10.3-dbl-cmplx_nomachopt
PETSC_DIR ?= /u/tbg/public/petsc
SLEPC_DIR ?= /u/tbg/public/slepc

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
 HDF5PATH = /afs/rzg/common/soft/hdf5/1.8.9/amd64_sles11/intel/12.1/impi/4.1.0/
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR) -Wl,-rpath,$(HDF5PATH)/lib
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(HAC),yes)
 #  HLST-ADIOS-CHECKPOINT interface
 ADIOS_DIR = /afs/ipp/itm/switm/lib/adios/adios-1.5.0-impi-O3/
 ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
 ADIOS_LIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
endif


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#e.g., specify list of files to be compiled w/o optimization in case of
#compiler bugs:
NOOPTLIST = 
