### architecture dependent GENE makefiles for C3SE beda cluster         ###
### http://www.c3se.chalmers.se/index.php/Hardware_Beda                 ###
###########################################################################
#									  #
# module load openmpi intel mkl						  #
#									  #
###########################################################################
### SWITCHES                                                            ###
###########################################################################

FFTLIB = mkl

#double precision should be default; single precision, however, might 
#be sufficient and much faster
PRECISION= double

#Switch on DEBUG mode in case you observed errors while running GENE
DEBUG= no

#Switch to yes if PETSC/SLEPC is installed (highly recommended!!)
SLEPC= no

#only required for the global code and full eigenvalue spectrum runs:
SCALAPACK = no

OPENMP = no
PERFLIB = none
FUTILS = no

#Provide an upper (RAM) memory limit for the code internal performance 
#optimization (w.r.t. MPI mapping & alternative routines)
PREPROC= -D'MB_PER_CORE=800'
#Note: the system itself often requires a significant fraction 
#      which should be taken into account

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS =

#FFT LIBRARY
ifeq ($(FFTLIB),mkl)
 INCPATHS += -I$(INCLUDE)
 MKLPATH=/beda/sw/apps/intel/Compiler/11.1/075/mkl
 MKLLIBPATH=$(MKLPATH)/lib/em64t
 LIBS += -L$(LIBRARY_PATH) -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a \
$(MKLLIBPATH)/libmkl_sequential.a $(MKLLIBPATH)/libmkl_core.a \
-Wl,--end-group -lpthread
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS +=
 LIBS +=
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
 PETSC_ARCH = DOUBLE_COMPLEX
 PETSC_DIR = /c3se/users/merz/software/petsc-3.0.0-p8
 SLEPC_DIR = /c3se/users/merz/software/slepc-3.0.0-p5
#slepc version < 3.0
# include $(SLEPC_DIR)/bmake/slepc_common
 include $(SLEPC_DIR)/conf/slepc_common
 INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
 LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
 PREPROC += -DWITHSLEPC
endif

#SCALAPACK is only necessary if you want to use the direct eigenvalue solver
ifeq ($(SCALAPACK),yes)
  LIBS +=
  PREPROC += -DWITHSCAL
endif

ifeq ($(FUTILS),yes)
 #  FUTILS and HDF5
 FUTILSDIR = $(EXTDIR)/futils/src
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
ifeq ($(DEBUG),yes)
  OPT = -O0
  FFLAGS += -C -traceback -debug inline-debug-info -warn all
  INCLUDE_SYMBOLS=yes
else
  OPT = -O2
#these switches would probably further improve performance;
#the intel compiler versions 12.0-12.1, however, do have some 
#bugs at the optimization levels
#  OPT= -O3 -ip -xSSSE3
endif

#FORTAN COMPILER FLAGS
FFLAGS = $(OPT) -module $(OBJDIR)

#LINKER (usually same as MPFC)
MPLD =$(MPFC)

#LDFLAGS
LDFLAGS = 

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################

#OPENMP
ifeq ($(OPENMP),yes)
FFLAGS +=
LDFLAGS += 
PREPROC +=-DWITHOMP
endif

MKFLAGS = $(FFLAGS)

ifeq ($(PRECISION),double)
FFLAGS += -r8 
PREPROC +=-DDOUBLE_PREC
endif

ifneq ($(PERFLIB),none)
 PREPROC += -DWITHPERF
endif

###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	export OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
	export MKL_SERIAL=yes;\
	mpiexec -n $(N_PES) ./$(EXEC)


###########################################################################
### Machine dependent compiling rules                                   ###
###########################################################################

ifeq ($(FFTLIB),mkl)
$(OBJDIR)/par_mod.o:           $(OBJDIR)/mkl_dfti.o
$(OBJDIR)/fourier_mkl.o:       $(OBJDIR)/mkl_dfti.o
$(OBJDIR)/mkl_dfti.o: 	       $(MKLPATH)/include/mkl_dfti.f90
	$(MPFC) $(MKLFLAGS) -c -o $@ $<  -module $(OBJDIR)
endif

#On Beda we need to use the mpimod wrapper:
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