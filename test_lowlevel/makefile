COMPILER = intel
BASEDIR = ..
OBJDIR = ./obj
SRCDIR = ../src

FFTLIB = fftw
USE_PERFLIB = perf
DEBUG = no
PRECISION=double
OPENMP = no
MKLVERSION=10.3
USE_CUDA_NONLIN=yes
INCLUDE_SYMBOLS=yes
COMPILER_REPORTS=no

MKLROOT=$(MKL_HOME)
OPTLEVEL=3

MACHMKFILE=

NVCC=nvcc -arch sm_20 -O3
COMP_CUDA_DEBUG = -g -G
CUDA_LIBS += -L$(CUDA_HOME)/lib64 -lcufft -lcudart

INCPATHS += -I$(OBJDIR) -I$(SRCDIR)
ifeq ($(FFTLIB),fftw)
	INCPATHS +=  -I$(FFTW_HOME)/include
ifeq ($(PRECISION),double)
	LIBS +=  -L$(FFTW_HOME)/lib -lfftw3
else
	LIBS +=  -L$(FFTW_HOME)/lib -lfftw3f
endif
endif

LIBS += $(MKL_SCALAPACK_LINKLINE)
ifeq ($(USE_CUDA_NONLIN),yes)
	LIBS += $(CUDA_LIBS)
endif

ifeq ($(PLATFORM),hpcff)
	CC=icc
	FC=mpif90
	COMP_MODDIR= -module $(OBJDIR)
	COMP_WARN=   -warn all
	COMP_MACHINE_OPT = -O4 -xHost
	COMP_DEBUG = -check all -traceback
	COMP_C_DEBUG = -traceback
	COMP_C_OPT = -O3 -xHost
	COMP_PRECISION = -real-size 64  -DDOUBLE_PREC
	COMP_WITHPERF = -DWITHPERF=1
	MKLLIBPATH = $(MKL_HOME)/lib/intel64
	COMP_LIBS = -Wl,--start-group -L$(MKLPATH) \
	-lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -Wl,--end-group

	FFLAGS = -m64
endif

ifeq ($(PLATFORM),power)
	FC=mpxlf95_r
	COMP_MODDIR=-qmoddir=$(OBJDIR)
	COMP_WARN = -qflag=I:E
	COMP_MACHINE_OPT = -O4 -qarch=pwr6 -qtune=pwr6
	FLTTRAP= -qflttrap=overflow:enable:invalid:zerodivide:NANQ
	COMP_DEBUG = -qcheck -qsigtrap -C -qinitauto=7FF7FFFF $(FLTTRAP)
	COMP_PRECISION = -qrealsize=8  -WF,-DDOUBLE_PREC
	COMP_LIBS = -L/usr/local/lib -lscalapack -llapack-essl -lblacsF77init -lblacs -lessl
	COMP_WITHPERF = -WF,-DWITHPERF

# Starting compiler flags
	FFLAGS = -q64
endif


#FORTAN COMPILER FLAGS
#FFLAGS += $(COMP_MODDIR) $(COMP_WARN) $(COMP_REPORT)
#LFLAGS += $(COMP_LFLAGS)
#LIBPATH =

ifeq ($(DEBUG),yes)
	FFLAGS += $(DEBUG_FFLAGS) $(DEBUG_NOUNUSED)
	CFLAGS += $(DEBUG_CFLAGS)
	CUDAFLAGS += $(COMP_CUDA_DEBUG)
else
ifeq ($(DEBUG),valgrind)
	FFLAGS += -O0
	CFLAGS += -O0
else
	FFLAGS += $(OPT_FFLAGS)
	CFLAGS += $(OPT_CFLAGS)
endif
endif

ifeq ($(COMPILER_REPORTS),yes)
	FFLAGS += $(REPORT_FFLAGS)
endif
ifeq ($(INCLUDE_SYMBOLS),yes)
	FFLAGS += -g
	CFLAGS += -g
endif

ifeq ($(OPENMP),yes)
	FFLAGS += -DWITHOMP -DWITHOMP_NONLIN $(OPENMP_FFLAGS)
	CFLAGS += $(OPENMP_CFLAGS)
endif

ifeq ($(USE_PERFLIB),hpm)
 HPM_DIR = $(IHPCT_BASE)
 INCPATHS += -I$(HPM_DIR)
 LIBS += -L$(HPM_DIR)/lib -L/afs/ipp/@sys/lib -lhpm_r -lperfhpm_r -lpmapi -lm -lxlsmp
endif
ifeq ($(USE_PERFLIB),perf)
 PREPROC += -DWITHPERF=2
 INCPATHS += -I$(PERFLIB_HOME)/include
 LIBS += -L$(PERFLIB_HOME)/lib -looperf -Wl,-rpath,$(PERFLIB_HOME)/lib
endif
ifeq ($(USE_PERFLIB),likwid)
 FFLAGS += -DWITH_LIKWID
 INCPATHS += -I$(LIKWID_HOME)/include
 LIBS += -L$(LIKWID_HOME)/lib -llikwid
endif
ifeq ($(USE_PERFLIB),papi)
 FFLAGS += -DWITHPAPI
 INCPATHS += -I$(PAPI_HOME)/include
 LIBS += -L$(PAPI_HOME)/lib -lpapi
endif
ifeq ($(USE_PERFLIB),none)
 OBJLIST += $(OBJDIR)/dummyperf.o
endif

ifeq ($(PRECISION),double)
	PREPROC += -DDOUBLE_PREC
	FFLAGS += $(DOUBLE_PRECISION_FFLAGS)
endif

FFLAGS += -heap-arrays 1024

FFLAGS += $(MODDIR_FFLAGS)
all:	checkpath nldrv

checkpath:
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)

include dependencies.mk
# dependencies on C files are done manually
$(OBJDIR)/df_nonlinear_cuda.o: $(OBJDIR)/cuda_nonlinearity.o
$(OBJDIR)/cuda_nonlinearity.o: 	$(OBJDIR)/cuda_overlap.o \
				$(OBJDIR)/cuda_fourier_cufft.o \
				$(SRCDIR)/cuda_kernels.cu \
				$(SRCDIR)/cuda_overlap.h
$(OBJDIR)/cuda_fourier_cufft.o:	$(SRCDIR)/cuda_overlap.h
$(OBJDIR)/cuda_overlap.o:	$(SRCDIR)/cuda_overlap.h

ifeq ($(USE_PERFLIB),likwid)
$(OBJDIR)/driver_nonlinear_term.o:	$(OBJDIR)/likwid.o
$(OBJDIR)/likwid.o: $(LIKWID_HOME)/include/likwid.f90
	$(FC) -c -o $@ $(MODDIR_FFLAGS) $<
endif
ifeq ($(USE_CUDA_NONLIN),yes)
LINK_OBJ += $(OBJDIR)/cuda_nonlinearity.o $(OBJDIR)/cuda_overlap.o $(OBJDIR)/cuda_fourier_cufft.o
endif

$(LINK_OBJ): makefile
NLDRVLIST = $(OBJDIR)/driver_nonlinear_term.o $(LINK_OBJ)

PMM_OBJLIST = perf_MatrixMatrixMultiplication.o BandedMatrix.o matrix.o Vector.o storefullmatrix.o processgrid.o \
	storebandedmatrix.o storevector.o ListObject.o

PMM_ASMLIST = $(subst .o,.s,$(PMM_OBJLIST))
PMM_SCLIST = $(subst .o,.sc,$(PMM_OBJLIST))

nldrv: $(NLDRVLIST)
	$(MPFC) -o nldrv $(FFLAGS) $(LFLAGS) $^ $(LIBS)

fixed:	small_matvec_fixed.o
	$(FC) -o fixed $(FFLAGS) $^ $(LIBS)

sc:	$(PMM_SCLIST)
bt:	blas_test.o
	$(FC) -o bt $(FFLAGS) $^ $(LIBS)

tcr:	test_cuda_routines.o
	$(NVCC) $(CUDAFLAGS) -o tcr $^ $(CUDA_LIBS)

tfc:	test_fourier_cufft.o fourier_cufft.o 
	$(CC) -o tfc $(CFLAGS) $^ -L$(FFTW_HOME)/lib -lfftw3 $(CUDA_LIBS)

pmm_asm: $(PMM_ASMLIST)

tfftwc: test_fftw_c.c
	$(CC) -o tfftwc -I$(FFTW_HOME)/include $(CFLAGS) $^ $(PERFLIB) -L$(FFTW_HOME)/lib -lfftw3

tfftw:	test_fftw.F90
	$(FC) -o tfftw -I$(FFTW_HOME)/include $(FFLAGS) $^ $(PERFLIB) -L$(FFTW_HOME)/lib -lfftw3 

derivmat_test:	grid1d.o processgrid.o \
		storefullmatrix.o matrix.o \
		derivative_matrix.o derivmat_test.o
	$(FC) -o derivmat_test  $(INCPATHS) $(FFLAGS) $^ $(LIBPATH) $(LIBS)

pmm:	$(PMM_OBJLIST)
	$(FC) -o pmm $^ $(LIBS)

tlo:	test_ListObject.o ListObject.o
	$(FC) -o tlo  $(INCPATHS) $(FFLAGS) $^ 

tbe:	test_boundary_exchange.o
	$(FC) -o tbe  $(INCPATHS) $(FFLAGS) $(TBE_OBJLIST)

tvm:	test_VectorModule.o
	$(FC) -o tvm  $(INCPATHS) $(FFLAGS) $(TVM_OBJLIST) $(LIBPATH) $(LIBS)

tmm:	test_MatrixModule.o 
	$(FC) -o tmm  $(INCPATHS) $(FFLAGS) $(TMM_OBJLIST) $(LIBPATH) $(LIBS)

tbmm:	storebandedmatrix.o test_StoreBandedMatrixModule.o
	$(FC) -o tbmm  $(INCPATHS) $(FFLAGS) $(TBMM_OBJLIST) $(LIBPATH) $(LIBS)

tbm:	BandedMatrix.o test_BandedMatrixModule.o
	$(FC) -o tbm $(INCPATHS) $(FFLAGS) $(TBM_OBJLIST) $(LIBPATH) $(LIBS)

mpi_test: mpi_send_recv_test.o
	$(FC) -o mpi_test $(FFLAGS) mpi_send_recv_test.o

fourier_cufft.o: fourier_cufft.c
	nvcc -c -o $@ $<

include ../makefiles/rules.mk

$(OBJDIR)/%.o: %.F90
	$(FC) $(FFLAGS) $(INCPATHS) $(PREPROC) -c -o $@ $<

#%.o:	%.F90
#	$(FC)  $(INCPATHS) $(FFLAGS) -c -o $@ $<

#%.o:	%.c
#	$(CC) $(INCPATHS) $(CFLAGS) -c -o $@ $<

#%.o:	%.cu
#	$(NVCC) $(CUDAFLAGS) -c -o $@ $<

#%.s:	%.F90
#	$(FC) $(INCPATHS) $(FFLAGS) -S -o $@ $<

#%.sc:	%.F90
#	$(FC) $(INCPATHS) $(FFLAGS) $<

.PHONY: run

# Dependencies
#$(OBJDIR)/driver_nonlinear_term.o: $(OBJDIR)/nonlinear_term.o \
				$(OBJDIR)/df_nonlinear_term.o

#nonlinear_term.o:	discretization.o \
			par_other.o \
			rhs_term.o
#par_in.o:		spectype.o


test_cuda_routines.o: cuda_kernels.cu reduction_kernel.cu


perf_MatrixMatrixMultiplication.o: matrix.o  Vector.o BandedMatrix.o

test_ListObject.o: ListObject.o

TMM_OBJLIST = test_MatrixModule.o matrix.o storefullmatrix.o processgrid.o \
	Vector.o Matrix_Vector.o storevector.o

test_MatrixModule.o:	matrix.o Vector.o Matrix_Vector.o
#Matrix_Vector.o:	matrix.o Vector.o

#matrix.o:	storefullmatrix.o
#storefullmatrix.o:	processgrid.o storevector.o

TVM_OBJLIST = test_VectorModule.o Vector.o processgrid.o storevector.o matrix.o storefullmatrix.o

test_VectorModule.o:	Vector.o
#Vector.o:		processgrid.o storevector.o matrix.o
#storevector.o:		processgrid.o

TBM_OBJLIST = test_BandedMatrixModule.o BandedMatrix.o processgrid.o storebandedmatrix.o matrix.o \
	Vector.o storevector.o storefullmatrix.o ListObject.o
test_BandedMatrixModule.o: BandedMatrix.o matrix.o Vector.o
#BandedMatrix.o:	processgrid.o storebandedmatrix.o matrix.o Vector.o

TBMM_OBJLIST = test_StoreBandedMatrixModule.o storebandedmatrix.o processgrid.o storefullmatrix.o ListObject.o storevector.o
#storebandedmatrix.o: processgrid.o storefullmatrix.o storevector.o ListObject.o

TBE_OBJLIST = test_boundary_exchange.o \
	boundary_exchange.o \
	setup_mpi_environment.o \
	dummyperf.o

test_boundary_exchange.o: boundary_exchange.o setup_mpi_environment.o BoundaryDescription.o

#boundary_exchange.o:    setup_mpi_environment.o \
			dummyperf.o \
			BoundaryDescription.o \
			boundary_exchange_x.o \
			boundary_exchange_z.o

#boundary_exchange_z.o:	boundary_exchange_x.o boundary_exchange_general.o box_data_module.o
#boundary_exchange_x.o:	boundary_exchange_general.o
#box_data_module.o:	grid1d.o par_mod.o
#par_mod.o:		spectype.o discretization.o coordinates.o par_other.o
#discretization.o:	par_in.o

mpi_io_test.o:	setup_mpi_environment.o



.PHONY:clean
clean:	
	-rm *.o 
	-rm $(OBJDIR)/*.o
	-rm $(OBJDIR)/*.mod
	-rm tbm
	-rm tfftw
	-rm tfftwc





