###########################################################################
### MAC OS / gfortran makefile                                          ###
###									###
### can also be used on debian systems:					###
### apt-get install gfortran openmpi-bin libopenmpi-dev 	        ###
### apt-get install libfftw3-3 liblapack-dev libatlas-base-dev          ###
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
SCALAPACK = no
USE_PERFLIB = none
DEBUG = no
WITH_CUTILS = no
###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR) 
INCPATHS += -I/usr/local/include -I/usr/include

#LIBRARIES AND LIBFLAGS
LIBS = -L/usr/local/gfortran/lib -L/usr/local/lib -L/usr/lib/ -lblas -llapack

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
 INCPATHS +=
 ifeq ($(PRECISION),double)
  LIBS += -lfftw3
 else
  LIBS += -lfftw3f
 endif
 LIBS +=
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# PETSC_ARCH = 
# PETSC_DIR = 
# SLEPC_DIR = 
 include $(SLEPC_DIR)/bmake/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS +=
  PREPROC += -DWITHSCAL
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif

###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpif90
ARCHIVE = ar r

FFLAGS =
ifeq ($(DEBUG),yes)
 OPT = -O0
 FFLAGS += -Waliasing -Wconversion -Wsurprising -Wunderflow -Wunused-parameter
 FFLAGS += -ffpe-trap=invalid,zero -fbacktrace
else
 OPT = -O3
endif
FFLAGS += $(OPT)

#FORTAN COMPILER FLAGS
FFLAGS += -J$(OBJDIR)  -fimplicit-none #-ffree-line-length-none

#PRECOMILER SWITCHES
PREPROC += -DMPI $(SVNDEF)

#LINKER (usually same as MPFC)
MPLD =$(MPFC)

#LDFLAGS
LDFLAGS = -lcblas

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#OPENMP
ifeq ($(OPENMP),yes)
FFLAGS +=
LDFLAGS += 
PREPROC +=-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qrealsize=8, etc.)
ifeq ($(PRECISION),double)
FFLAGS += -fdefault-real-8
PREPROC +=-DDOUBLE_PREC
endif


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS);
	MALLOC_CHECK_=0;\
	mpirun -np $(N_PES) ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#gfortran -O3 optimization causes crashes in standard test 3,
#perf_vec(4)=2
$(OBJDIR)/collisions.o:	$(SRCDIR)/collisions.F90
			@echo "Compiling " $<
			$(MPFC) $(FFLAGS) $(INCPATHS) $(PREPROC) -c -O2 -o $@ $<

#Do we need to use the mpimod wrapper on Macs?
MPIDEPOBJ = \
	$(OBJDIR)/BoundaryDescription.o\
	$(OBJDIR)/boundary_exchange_general.F90\
	$(OBJDIR)/boundary.o\
	$(OBJDIR)/checkpoint.o\
	$(OBJDIR)/comm.o\
	$(OBJDIR)/field_solve_xky.o\
	$(OBJDIR)/gene.o\
	$(OBJDIR)/gyro_average_xky.o\
	$(OBJDIR)/storebandedmatrix.o\
	$(OBJDIR)/storefullmatrix.o\
	$(OBJDIR)/storevector.o

$(MPIDEPOBJ): $(OBJDIR)/mpimod.o
