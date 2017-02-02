###########################################################################
### architecture dependent GENE makefiles, palma cluster version        ###
###########################################################################
#									  #
# the following modules have to be loaded first:		          #
# module load fftw							  #
# module load mkl							  #
# module load petsc							  #
# module load slepc							  #
#								          #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

#COMPILER: INTEL
COMPTYPE = INTEL

# FFTLIB - needed in /src/files.mk
# set to: mkl or fftw
FFTLIB = fftw
PRECISION= double
OPENMP = no

SLEPC = yes
# has to be switched on
SCALAPACK = no
DEBUG = no
INCLUDE_SYMBOLS = no
COMPILER_REPORTS = no
USE_PERFLIB = none
FUTILS = no
PRODRUN = yes
WITH_CUTILS = yes

#memory per core
PREPROC = -D'MB_PER_CORE=1801'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

MKLLIBPATH = $(MKLPATH)/lib/em64t

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS

### FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 MKLFLAGS += 
 INCPATHS += -I$(MKLPATH)/include/em64t -L$(MKLPATH)
# LIBS += -lmkl_em64t -Bstatic -lguide -Bdynamic -lpthread
endif
ifeq ($(FFTLIB),fftw)
  ifeq ($(PRECISION),double)
   INCPATHS += -I$(FFTW_INCLUDE_DIR)
  LIBS += -L$(FFTW_LIB_DIR) -lfftw3
  else
  INCPATHS += -I$(FFTW_INCLUDE_DIR)
  LIBS += -L$(FFTW_LIB_DIR) -lfftw3f
  endif
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
# PETSC_ARCH = linux-gnu-c-opt
 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes) 
  LIBS += -Wl,--start-group $(MKLLIBPATH)/libmkl_scalapack_lp64.a \
 $(MKLLIBPATH)/libmkl_blacs_intelmpi_lp64.a $(MKLLIBPATH)/libmkl_intel_lp64.a \
 $(MKLLIBPATH)/libmkl_sequential.a $(MKLLIBPATH)/libmkl_core.a -Wl,--end-group \
 -lpthread 
  PREPROC += -DWITHSCAL
else
##Insert BLAS library
 MKLLIBPATH = $(MKLPATH)/lib/em64t
 LIBS += -L$(MKLLIBPATH)
 LIBS += -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a \
 $(MKLLIBPATH)/libmkl_sequential.a $(MKLLIBPATH)/libmkl_core.a \
 -Wl,--end-group -lpthread
endif

#ifeq ($(FUTILS),yes)
# #  FUTILS and HDF5
# FUTILSDIR = $(EXTDIR)/$(MACHINE)/futils/src
# HDF5PATH = /afs/ipp/home/t/tbg/soft/@sys/hdf5
# #HDF5PATH = $(HDF5_HOME)
# HDF5_LIBPATH  = -L$(HDF5PATH)/lib -L$(FUTILSDIR)
# HDF5_LIBS = -lfutils -lhdf5_fortran -lhdf5 -lz 
#
# LIBS += $(HDF5_LIBPATH) $(HDF5_LIBS)
# PREPROC +=  -DWITHFUTILS
# INCPATHS += -I$(FUTILSDIR)
#
#endif


ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif

#ifeq ($(USE_PERFLIB),perf)
# LIBS += -lperf_r
#endif

###########################################################################
### COMPILER & COMPILER FLAGS                                           ###
###########################################################################

#PRECOMILER SWITCHES
PREPROC += $(F2003_MISSING) -DLINUX -DMPI $(SVNDEF)
CPREPROC += -DLINUX

MPFC = mpiifort
ARCHIVE = ar r

FFLAGS =

ifeq ($(COMPTYPE),INTEL)
 ifeq ($(DEBUG),yes)
  OPT = -O0
  FFLAGS += -C -traceback -debug inline-debug-info -warn all
#'-check noarg_temp_created' suppresses warnings about temporary arrays during runtime 
  INCLUDE_SYMBOLS=yes
 else
  ifeq ($(PRODRUN),yes)
   OPT= -O3 -xO -us -pad -funroll-loops -ip -ipo
  else
   OPT = -O3 
  endif
  OPT += -vec_report0
 endif

 ifeq ($(INCLUDE_SYMBOLS),yes)
   FFLAGS += -g
 endif

 ifeq ($(COMPILER_REPORTS),yes)
   FFLAGS += -opt-report -opt-report-phase=all \
	     -opt-report-file=gene_report.opt \
	     -vec-report3
 endif
 MKLFLAGS = -module $(OBJDIR)

 MPLD =$(MPFC)
 LDFLAGS = -lsvml
endif

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
#OPENMP
ifeq ($(OPENMP),yes)
 ifeq ($(COMPTYPE),INTEL)
  MKLFLAGS += -openmp -fpp
  LDFLAGS += -openmp
 endif
PREPROC +=-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qautodbl, etc.)
ifeq ($(PRECISION),double)
 ifeq ($(COMPTYPE),INTEL)
  FFLAGS += -r8
 endif
PREPROC +=-DDOUBLE_PREC
endif

FFLAGS += $(MKLFLAGS) $(OPT)


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun 
run:	$(EXECDIR)/$(EXEC)
	cd $(RUNDIR);\
#	LD_LIBRARY_PATH=\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS)\
#	MKL_SERIAL=YES\
	mpirun -np $(N_PES) ./$(EXEC)

###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

ifeq ($(FFTLIB),mkl)
 $(OBJDIR)/par_mod.o:		$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/fourier_mkl.o:	$(OBJDIR)/mkl_dfti.o
 $(OBJDIR)/mkl_dfti.o: $(MKLPATH)/include/mkl_dfti.f90
	$(MPFC) $(MKLFLAGS) -c -o $@ $<
endif

NOOPTLIST= 

NOOPTFULLOBJ = $(addprefix $(OBJDIR)/,$(NOOPTLIST))

$(NOOPTFULLOBJ): $(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(MPFC) $(NOOPTFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<
