###########################################################################
### architecture dependent GENE makefiles, gp cluster version          ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load impi/4.0.3 gcc/4.6.3 mkl/10.3 fftw petsc slepc				  # 
# module load hdf5-mpi perflib/2.0                                               #
#								          #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

COMPTYPE = GNU

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = mkl
PRECISION= double
OPENMP = no

SLEPC = no
# has to be switched on
SCALAPACK = yes
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS=no
USE_PERFLIB = none
FUTILS = no
PRODRUN = no
MKLVERSION = 10.2
USE_MEMORY_CHECKER=none

#memory per core
PREPROC = -D'MB_PER_CORE=1900'
FCK_PREPROC = 'MP_PER_CORE=1900'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#MKLPATH = $(MKL_HOME)
#MKLLIBPATH = $(MKL_HOME)/lib/em64t
MKLROOT=$(MKLPATH)/../..

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

FCK_INCPATHS = $(OBJDIR):.:$(SRCDIR)
#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 MKLFLAGS += 
 INCPATHS +=
endif
ifeq ($(FFTLIB),fftw)
#FFTW seems to be faster than MKL when single precision is used
#further tests needed
   INCPATHS += -I$(FFTW_HOME)/include
   FCK_INCPATHS := $(FCK_INCPATHS):$(FFTW_HOME)/include
 ifeq ($(PRECISION),double)
   LIBS += -L$(FFTW_HOME)/lib -lfftw3
 else
   LIBS += -L$(FFTW_HOME)/lib -lfftw3f
 endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  ifeq ($(PRECISION),double)
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(SLEPC_LIB)
   PREPROC += -DWITHSLEPC
  else
   SLEPC=no
  endif
endif

ifeq ($(SCALAPACK),yes)
ifeq ($(MKLVERSION),10.3)
# Intel MKL 10.3 sequential
	LIBS += -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 \
		-lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
		-lmkl_blacs_intelmpi_lp64 -lpthread -lm
	MKLFLAGS =  -m64 
	INCPATHS += -I$(MKLROOT)/include
	LDFLAGS += -Wl,-rpath,$(MKLROOT)/lib/intel64
else
# Intel MKL 10.2 sequential
	LIBS += $(MKLROOT)/lib/em64t/libmkl_scalapack_lp64.a \
		$(MKLROOT)/lib/em64t/libmkl_solver_lp64_sequential.a \
		-Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_gf_lp64.a \
		$(MKLROOT)/lib/em64t/libmkl_sequential.a \
		$(MKLROOT)/lib/em64t/libmkl_core.a \
		$(MKLROOT)/lib/em64t/libmkl_blacs_intelmpi_lp64.a \
		-Wl,--end-group -lpthread -lm

	MKLFLAGS += -m64 
	INCPATHS += -I$(MKLROOT)/include
	#LDFLAGS += -Wl,-rpath,$(MKLROOOT)/lib/em64t
endif
  PREPROC += -DWITHSCAL
  FCK_PREPROC := $(FCK_PREPROC):WITHSCAL
else
#Insert BLAS library
LIBS += -L$(MKLROOT)/lib/intel64 \
	-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
MKLFLAGS += -m64 
INCPATHS += -I$(MKLROOT)/include
LDFLAGS += -Wl,-rpath,$(MKLROOT)/lib/intel64
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

MEMCHECK_OBJ=
ifeq ($(USE_MEMORY_CHECKER),own)
	MEMCHECK_OBJ = $(OBJDIR)/set_malloc_hooks.o
	PREPROC += -DMEMORY_CHECKER
endif

###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################
#PRECOMILER SWITCHES
PREPROC += $(F2003_MISSING) -DLINUX -DMPI
FCK_PREPROC := $(FCK_PREPROC):LINUX:MPI
CPREPROC += -DLINUX

FFLAGS =

# first set default for MPFC
MPFC = mpif90
ifeq ($(USE_PERFLIB),scalasca)
	MPFC := scalasca -instrument -user -f90 -rcfile ./opari.rc -- $(MPFC)
	INCPATHS += -I/afs/ipp/home/t/tisd/gp_inst/include
	FFLAGS += -fpp -Wp,-P
endif
ifeq ($(USE_PERFLIB),ompp)
# the preprocessor must not write #line directives
# and also the quotes in the SVN statement is switched off
	MPFC := kinst-ompp $(MPFC)
	FFLAGS += -fpp -Wp,-P
else
#	PREPROC += $(SVNDEF)
endif

ifeq ($(USE_PERFLIB),perf)
 LIBS += $(PERFLIB_TS)
 PREPROC += -DWITHPERF=1
endif

ifeq ($(USE_PERFLIB),vtrace)
 FFLAGS += -tcollect
 LDFLAGS += -tcollect
# PREPROC += -DVTRACE
 INCPATHS += -I$(ITAC_HOME)/include
endif

ARCHIVE = ar r
MPLD = $(MPFC)

ifeq ($(COMPTYPE),GNU)
  ifeq ($(DEBUG),yes)
    FFLAGS += -Wall -Warray-temporaries -ffpe-trap=invalid,zero,overflow,underflow
    FFLAGS += -fcheck=all
  else
    FFLAGS += -O3 -march=corei7
  endif

  ifeq ($(INCLUDE_SYMBOLS),yes)
    FFLAGS += -g -fbacktrace
  endif
endif


#LIBS += -L$(GCC_HOME)/lib64 -lgfortran


###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#OPENMP
ifeq ($(OPENMP),yes)
 ifeq ($(COMPTYPE),INTEL)
  MKLFLAGS += -openmp 
  LDFLAGS += -openmp
 endif
 ifeq ($(COMPTYPE),LAHEY)
  MKLFLAGS += --openmp
  LDFLAGS += --openmp
 endif 
PREPROC +=-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qautodbl, etc.)
ifeq ($(PRECISION),double)
 ifeq ($(COMPTYPE),GNU)
  FFLAGS += -fdefault-real-8 -fdefault-double-8
 endif
PREPROC +=-DDOUBLE_PREC
FCK_PREPROC := $(FCK_PREPROC):DOUBLE_PREC
endif

MODDIR = -J$(OBJDIR)
FFLAGS += -std=f2003 -fall-intrinsics $(MODDIR) $(MKLFLAGS) $(OPT)
#FFLAGS += $(MKLFLAGS) $(OPT)


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun 
run:	$(EXECDIR)/$(EXEC)
	cd $(RUNDIR);\
	export OMP_NUM_THREADS=1;\
	export MKL_SERIAL=yes;\
	mpiexec -n $(N_PES) -env MV2_USE_SRQ 0 ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

ifeq ($(FFTLIB),mkl)
 $(OBJDIR)/par_mod.o:		$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/fourier_mkl.o:	$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/mkl_dfti.o: $(MKLROOT)/include/mkl_dfti.f90
	$(MPFC) $(MKLFLAGS) -c -o $@ $<
endif

NOOPTLIST= 

NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<


#On BOB we need to use the mpimod wrapper:
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

forcheck:
	forchk -define $(FCK_PREPROC) -I $(FCK_INCPATHS) -l mylistfile -ff -i4 -dp -allc -decl -nshsrc -ancmpl -anref -shmod $(F90FULLSRC) -lib $(FCKDIR)/share/forcheck/MPI.flb
