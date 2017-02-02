###########################################################################
### architecture dependent GENE makefile, ANSELM version                ###
### https://docs.it4i.cz/anselm-cluster-documentation                   ###
###########################################################################

###########################################################################
### List required environment settings/modules here:                    ###
### module load PrgEnv-intel fftw3                                      ###
###########################################################################

###########################################################################
### Note: In the following, space characters must not be added at the   ###
###       end of variable names!                                        ###
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = intel
#COMPILER = $(shell echo $(PE_ENV) | tr A-Z a-z)

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP = SandyBridge

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = mpirun -np $(N_PES) -hostfile $(PBS_NODEFILE) ./$(EXEC)

###########################################################################
### SWITCHES                                                            ###
###########################################################################


# FFTLIB - needed in /src/files.mk and /makefiles/rules.mk
# set to: mkl, fftw or essl
FFTLIB = fftw

#double precision should be default; single precision, however, might 
#be sufficient and much faster
PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= yes

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = yes

#OPENMP might be important in future GENE releases again
#Currently, pure MPI is most likely the best choice
OPENMP = no

#performance profiling library 
# * switch to none for normal simulations (!)
# * possible choice distributed with GENE is USE_PERFLIB=FR
USE_PERFLIB = none

# FUTILS and HDF5 are required, e.g. for some geometry interfaces
FUTILS = no

#Provide an upper (RAM) memory limit in MB per core for the code internal 
#performance optimization (w.r.t. MPI mapping & alternative routines)
MB_PER_CORE=3500
#Note: the system itself often requires a significant fraction 
#      which should be taken into account

#include symbols or compiler reports for debugging
#(will automatically be activated by DEBUG<>no)
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


MKLVERSION=13.5


#FFT LIBRARY
ifeq ($(FFTLIB),mkl)
#set MKLINCLUDE PATH here
 MKLINCLUDE = -I$(MKLROOT)/include
 INCPATHS += -I$(MKL_INC_DIR)
 LIBS += -Wl,--start-group $(MKLROOT)/lib/$(ARCH)/libmkl_intel
_lp64.a\
      $(MKLROOT)/lib/$(ARCH)/libmkl_core.a\
      $(MKLROOT)/lib/$(ARCH)/libmkl_sequential.a\
      -Wl,--end-group -lpthread -lm
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(FFTW3_INC_DIR)
 LIBS += -L$(FFTW3_LIB_DIR) -lfftw3 -lfftw3f
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH = icc-impi-mkl-opt
 PETSC_DIR = /home/pr1g0010/soft/petsc/
 SLEPC_DIR = /home/pr1g0010/soft/slepc/

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
else
  LIBS += $(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_INC_DIR)/../
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
