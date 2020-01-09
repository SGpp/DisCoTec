###########################################################################
### architecture dependent GENE makefiles, Ubuntu+gfortran laptop version #
###########################################################################
#									  #
# the following packages have to be installed first:		          #
# apt-get install gfortran openmpi-bin libopenmpi-dev 	        	  #
# apt-get install libfftw3-3 liblapack-dev libatlas-base-dev		  #
#								          #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPILER=gnu
CHIP=

MPRUN = export OMP_NUMTHREADS=1; export MKL_SERIAL=yes;\
      mpiexec -np $(N_PES) ./$(EXEC)

FFTLIB = fftw
PRECISION= double
OPENMP = no

MPI_IO = no
SLEPC = yes 
SCALAPACK = no
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = no

#memory per core
PREPROC = -D'MB_PER_CORE=env'
FCK_PREPROC = 'MP_PER_CORE=1000'
COMPILE_MPI_MOD=yes

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################


#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)
INCPATHS += -I/usr/local/include -I/usr/include
SLEPC_DIR=/home/sccs/personal/downloads/slepc-3.5.3
PETSC_DIR=/home/sccs/personal/downloads/petsc-3.5.3
PETSC_ARCH=arch-linux-gnu-c-debug-complex


SLEPC_INCLUDE=$(SLEPC_DIR)/$(PETSC_ARCH)/include
PETSC_INCLUDE=$(PETSC_DIR)/$(PETSC_ARCH)/include
#LIBRARIES AND LIBFLAGS
LIBS = -L/usr/local/gfortran/lib -L/usr/local/lib -L/usr/lib/ -L$(PETSC_DIR)/$(PETSC_ARCH)/lib/ -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib/ -llapack -lblas -lslepc -lpetsc

### FFT LIBRARY
ifeq ($(FFTLIB),fftw)
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
 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_INCLUDE) $(SLEPC_INCLUDE)
 LIBS += -lslepc -lpetsc
 PREPROC += -DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
  LIBS += $(MKL_SCALAPACK_LINKLINE)
  PREPROC += -DWITHSCAL
  FCK_PREPROC := $(FCK_PREPROC):WITHSCAL
else
#Insert only BLAS library
  LIBS += #$(MKL_BLAS_LINKLINE)
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = $(HDF5_HOME)
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)
endif

ifeq ($(USE_PERFLIB),perf)
 LIBS += $(PERFLIB_TS)
 PREPROC += -DWITHPERF=1
else
ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif
endif

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

ifeq ($(FFTLIB),mkl)
$(OBJDIR)/mkl_dfti.o: $(MKLROOT)/include/mkl_dfti.f90
	$(MPFC) $(MODDIR_FFLAGS) -c -o $@ $<
endif

NOOPTLIST= 

NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<
