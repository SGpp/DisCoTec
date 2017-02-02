###########################################################################
### template for architecture dependent GENE makefiles                  ###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = essl
#FFTLIB = fftw
PRECISION= double
OPENMP = no

SLEPC= no
# perflib can take values of hpm or none
USE_PERFLIB = none
BITMODE = 64

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
#LIBS = -L/gpfs/apps/BLAS/blas-1.0.0/64/lib -lblas

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS += -L/usr/lib64 -lessl
 PREPROC += -WF,-DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I/gpfs/apps/FFTW/64/include
 LIBS += -L/gpfs/apps/FFTW/64/lib
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
 PETSC_DIR = /gpfs/apps/PETSC/2.3.3
 SLEPC_DIR = /gpfs/apps/SLEPC/SRC/slepc-2.3.3
 include  $(SLEPC_DIR)/bmake/slepc_common
 INCPATHS += $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS += $(PETSC_LIB) $(SLEPC_LIB) 
 PREPROC += -WF,-DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
  LIBS +=
  PREPROC += -WF,-DWITHSCAL
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -WF,-DWITHPERF
endif
ifeq ($(USE_PERFLIB),hpm)
 LIBS += -L/gpfs/apps/IHPCT/2.1/lib64 -lhpm -llicense -lshpm
endif

###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpif90 

#ARCHIVE command
ARCHIVE = ar -X$(BITMODE) r

#C COMPILER
CC = xlc 

#FORTAN COMPILER FLAGS
FFLAGS = -q$(BITMODE) -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto \
-qmoddir=$(OBJDIR) -qflag=I:I

#PRECOMILER SWITCHES
PREPROC += -WF,-DAIX

#LINKER (usually same as MPFC)
MPLD =$(MPFC)

#LDFLAGS
LDFLAGS = -q64 

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#OPENMP
ifeq ($(OPENMP),yes)
FFLAGS += -qsmp=omp
LDFLAGS += -qsmp=omp
PREPROC +=-WF,-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qrealsize=8, etc.)
ifeq ($(PRECISION),double)
FFLAGS +=-qrealsize=8
PREPROC +=-WF,-DDOUBLE_PREC
endif


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)


#------- mpirun command (needed in testsuite, scanscript)
run:	$(EXECDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS)\
	srun ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

