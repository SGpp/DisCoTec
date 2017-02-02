###########################################################################
### AKKA version of machine dependent GENE makefiles                    ###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw
PRECISION = double
OPENMP = no

SLEPC = no
USE_PERFLIB = none

#memory per core
#please fill and uncomment the following line
#PREPROC= -D'MB_PER_CORE='

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS += -lgoto -llapack -lgoto
LDFLAGS += 


#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
#this is FFTW 3.2.1, configured with psc and openmpi/psc modules loaded and
#./configure CC=mpicc --enable-sse2 --prefix=/home/f/flmerz/software/fftdbl_lib
#can be removed once the fftw module is recompiled
 MYFFTW=/home/f/flmerz/Public/fftdbl_lib
 INCPATHS += -I$(MYFFTW)/include
 LIBS += -L$(MYFFTW)/lib -lfftw3

# this can be used again once the module is fixed
# INCPATHS += $(FFTW_INCLUDE)
# LIBS += -lfftw3 -lfftw3f $(FFTW_LDFLAGS)  
# LDFLAGS += $(FFTW_LDFLAGS)
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH = linux-psc-c-opt
 PETSC_DIR = /home/s/skymandr/Public/software/petsc/3.0.0-p5
 SLEPC_DIR = /home/s/skymandr/Public/software/slepc/3.0.0-p3 
 include $(SLEPC_DIR)/bmake/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif

###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpif90

#ARCHIVER
ARCHIVE = ar r

#FORTAN COMPILER FLAGS
# FFLAGS = 
FFLAGS = -O3 -module $(OBJDIR)

#PRECOMILER SWITCHES
PREPROC += -DMPI $(SVNDEF)

#LINKER (usually same as MPFC)
MPLD = $(MPFC)

#LDFLAGS
LDFLAGS = -L/usr/local/lib

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#OPENMP
ifeq ($(OPENMP),yes)
# OMP_NUM_THREADS = 1
 LDFLAGS += -openmp
 FFLAGS += -openmp -fpp
 PREPROC +=-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qrealsize=8, etc.)
ifeq ($(PRECISION),double)
 FFLAGS += -r8
 PREPROC += -DDOUBLE_PREC
endif


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
	mpiexec -np $(N_PES) ./$(EXEC)

###########################################################################
### Machine dependent ompiling rules                                    ###
###########################################################################

