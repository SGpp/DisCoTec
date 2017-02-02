###########################################################################
### PPPL cluster              		     				###
### researchcomputing.pppl.gov                                          ###
###########################################################################

###########################################################################
#                                                                         #
# Last Tested 5/21/2014 S. Lazerson (lazerson@pppl.gov)                   #
# load the following modules before compilation with *pgi*                #
# module purge                                                            #
# module load ppplcluster                                                 #
# module load pgf/13.6                                                    #
# module load openmpi/1.4.3-qlogic acml superlu                           #
# module load blacs metis parmetis superlu_dist                           #
# module load fftw scalapack                                              #
# module load slepc_complex petsc_complex                                 #
# module load subversion                                                  #
#                                                                         #
# for launcher: module load python/3.2.2                                  #
# for post-processing: module load gnuplot idl                            #
###########################################################################

###########################################################################
### BASIC settings                                                      ###
###########################################################################

#set compiler (check /makefiles/compilers/*.def for possible choices)
COMPILER = pgi

MPFC ?= mpif90
MPCC ?= mpicc

FC ?= pgf90
CC ?= pgcc

GENE_COMPILE_VERBOSE = 0

#set chip for proper chip optimization flags (optional)
#check possible choices in /makefiles/compilers/$COMPILER.def
CHIP = 

#ARCHIVE command
ARCHIVE = ar r

#MPI command to run the executable $(EXEC) with $N_PES MPI processes
MPRUN = mpiexec -n $(N_PES) ./$(EXEC)

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
MB_PER_CORE=1750
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
INCPATHS = -I$(OBJDIR) -I$(SRCDIR) 

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LDFLAGS = 
LIBS = -L$(ACMLHOME)/pgi64/lib/ -lacml

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
 INCPATHS += -I$(FFTWHOME)/include
 LIBS += -L$(FFTWHOME)/lib -lfftw3
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  include  $(SLEPC_DIR)/conf/slepc_common
  INCPATHS += $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
  LIBS     += $(SLEPC_LIB) 
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  INCPATHS += -I${BLACS_HOME}/include
  LIBS += -L${SCALAPACK_HOME} -lscalapack 
  LIBS += -L${BLACS_HOME}/lib -lmpiblacs -lmpiblacsCinit -lmpiblacsF77init
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
