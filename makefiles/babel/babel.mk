###########################################################################
### architecture dependent GENE makefile for RZG BlueGene/P             ###
###########################################################################

###########################################################################
### SWITCHES                                                            ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw or essl
FFTLIB = essl
PRECISION = double
OPENMP = no

SLEPC= yes
# USE_PERFLIB set to perf, none, hpm or scalasca
USE_PERFLIB = none
MASSLIB = no
DEBUG = no
PRODRUN = yes
SCALAPACK = yes

#memory per core
PREPROC= -WF,-D'MB_PER_CORE=350'

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
#Add include paths for BLAS routines (can be LAPACK,ESSL,MKL,etc.)
BGP_SYS = /bgsys/drivers/ppcfloor/comm
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR) -I$(BGP_SYS)/include

#LIBRARIES AND LIBFLAGS
#Insert BLAS library
LIBS = -L$(BGP_SYS)/lib

#FFT LIBRARY
#specify at least one of the following choices
ifeq ($(FFTLIB),mkl)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS += -L/opt/ibmmath/lib -lesslbg
 PREPROC += -WF,-DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS +=
 LIBS +=
endif

###########################################################################
# ADDITIONAL LIBRARIES (set via SWITCH in header)			  #
###########################################################################

ifeq ($(SLEPC),yes)
  PETSC_DIR=/workgpfs/deisa/rzg33gen/rzg33fsj/PETSc
  SLEPC_DIR=/workgpfs/deisa/rzg33gen/rzg33fsj/SLEPc
  PETSC_ARCH =linux-gnu-c-opt
  include $(SLEPC_DIR)/conf/slepc_common
  INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
  #LIBS +=$(PETSC_LIB) $(SLEPC_LIB)
  LIBS += $(SLEPC_LIB)
  PREPROC += -WF,-DWITHSLEPC
endif

ifeq ($(SCALAPACK),yes)
  LIBS += -L/usr/local/lib -lscalapack -llapack -lblacsF77init -lblacs 
  LIBS += -L/opt/ibmmath/lib -lesslbg
  PREPROC += -WF,-DWITHSCAL
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -WF,-DWITHPERF
endif

ifeq ($(USE_PERFLIB),perf)
 LIBS += -L/usr/local/lib
endif

ifeq ($(MASSLIB),mass)
 MASS_LIB = -lmass
endif

###########################################################################
### COMPILER & COMPILER FLAGS       				  	###
###########################################################################

#FORTRAN COMPILER
ifeq ($(USE_PERFLIB),scalasca)
  MPFC = scalasca -instrument -user mpixlf90_r
else
  MPFC = mpixlf90_r
endif
ARCHIVE = ar r 

#FORTAN COMPILER FLAGS
ifeq ($(PRODRUN),yes)
 OPTLEVEL = 4
 IPA= -qipa
else
 OPTLEVEL = 3
endif
ifeq ($(DEBUG),yes)
 OPT= -qcheck -g -O0 -C
else
 OPT= -O$(OPTLEVEL) $(IPA) -Q -qarch=450d -qtune=450\
	-qunroll -qalias=:nopteovrlp:nointptr:std
endif
FFLAGS = $(OPT) -qmaxmem=-1 -qmoddir=$(OBJDIR) -qflag=I:I

#PRECOMILER SWITCHES
PREPROC += -WF,-DAIX
PREPROC += -WF,$(SVNDEF)

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
PREPROC += -WF,-DWITHOMP
endif

#Specify compiler flag for conversion to double precision
#(e.g. -r8, -qrealsize=8, etc.)
ifeq ($(PRECISION),double)
FFLAGS += -qautodbl=dbl4
PREPROC += -WF,-DDOUBLE_PREC
endif


###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS)\
	mpirun -np $(N_PES) -mode VN -exe ./$(EXEC)

###########################################################################
### Machine dependent compiling rules                                    ###
###########################################################################

#On Babel we need to use the mpimod wrapper:
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


# The following compiling rule is due to a compiler bug. If used -O3 
# the compiler does make an error so we need -qstrict also.
# This rule can be deleted, if the compiler has been repaired.
$(OBJDIR)/vel_space.o: $(SRCDIR)/vel_space.F90
	@echo "Special compiling rule for " $<
	$(MPFC) $(FFLAGS) -qstrict $(INCPATHS) $(PREPROC) -c -o $@ $<

