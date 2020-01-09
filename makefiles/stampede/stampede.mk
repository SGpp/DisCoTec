###########################################################################
### template for architecture dependent GENE makefiles                  ###
###########################################################################
#module load fftw3
#export MACHINE=stampede
###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = fftw
PRECISION= double
OPENMP = no

SLEPC= yes
SCALAPACK = yes
WITHNAG = no
USE_PERFLIB = none

#Switch GENE diagnostics output to MPI_IO:
MPI_IO= no

# FUTILS and HDF5 are required for some geometry interfaces
FUTILS = no

#insert memory per core and uncomment the following line
PREPROC= -D'MB_PER_CORE=1500'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR) -I$(MKLROOT)/include  

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 MKLFLAGS +=
 INCPATHS += 
 LIBS += 
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS += -I$(TACC_FFTW3_INC)
 ifeq ($(PRECISION),double)
  LIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpi -lfftw3
 else
  LIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpif -lfftw3f
 endif


endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(WITHNAG),yes)
 INCPATHS +=
 LIBS +=
#DEFAULT: use NAG-FORTRAN 90
 PREPROC +=-DWITHNAGF90
#Alternative: use NAG_FORTRAN-F77
#PREPROC +=-DWITHNAGF77 
endif

ifeq ($(SLEPC),yes)

 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS += $(SLEPC_INCLUDE)
 INCPATHS += -I${PETSC_DIR}/include -I${PETSC_DIR} -I${PETSC_DIR}/${PETSC_ARCH}/include
 LIBS += $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC

endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  INCPATHS += -I$(MKLROOT)/include  
  LIBS += $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm
  PREPROC += -DWITHSCAL
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
 HDF5PATH = 
 HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
 HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 

 LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
 PREPROC +=  -DWITHFUTILS
 INCPATHS += -I$(FUTILSDIR)

endif


###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
MPFC = mpif90

#ARCHIVE command
#Note: On IBM machines you might need to add -X$(BITMODE) to the "ar" call
ARCHIVE = ar r

#FORTAN COMPILER FLAGS
FFLAGS = -xhost

#PRECOMILER SWITCHES
PREPROC += -DMPI $(SVNDEF)

ifeq ($(DEBUG),yes)
  OPT = -O0
  FFLAGS += -g 
else
#Should test others as well
  OPT = -O2 
endif
FFLAGS += $(OPT) 

MKLFLAGS = -module $(OBJDIR)

#LINKER (usually same as MPFC)
LD =$(MPFC)

#LDFLAGS
LDFLAGS =

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#MPI-IO (only with full MPI2 library)
ifeq ($(MPI_IO),yes)
 PREPROC += -DWITHMPI_IO
endif

#OPENMP
ifeq ($(OPENMP),yes)
FFLAGS += -mp
LDFLAGS +=-mp
PREPROC +=-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qrealsize=8, etc.)
ifeq ($(PRECISION),double)
FFLAGS += -r8
PREPROC += -DDOUBLE_PREC
#PREPROC += -D'erf=derf'
endif

FFLAGS += $(MKLFLAGS)

###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(RUNDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS)\
	ibrun ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

#>>>> Try this:
ifeq ($(FFTLIB),mkl)
 $(OBJDIR)/par_mod.o:		$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/fourier_mkl.o:	$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/mkl_dfti.o: $(TACC_MKL_DIR)/include/mkl_dfti.f90
	$(MPFC) $(MKLFLAGS) -c -o $@ $<
endif

