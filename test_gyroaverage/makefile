COMPILER = intel
OBJDIR = ./obj_$(MACHINE)
SRCDIR = ../src
CFLAGS =
FFLAGS = 
FFTLIB = fftw
USE_PERFLIB = none
DEBUG = yes
PRECISION=double
MKLVERSION=10.3
MKLROOT=$(MKL_HOME)

include ../makefiles/compilers/$(COMPILER).def

#this default target has to be called first
# it runs the preprocessor with the variables set in the
# $(MACHINE).mk file, which is local to the test_gyroaverage directory
# then it runs the get_dependencies script, which builds the dependencies
# of the fortran modules. They are written into dependencies.mk.
# As last part this makefile is again called (but now the dependencies are right)
# and the all target is used.

# The gene mastermakefile and the makefile structure of GENE is
# only used for the preprocessing. The actual compilation of the
# test programs is steered by this makefile here.
preprocess:
	touch make_local
	make -f ../makefile preproc
	../tools/get_dependencies.py -r test_gyroaverage_df preproc/* test_gyroaverage.F90 test_pmm.F90 test_gyroaverage_io.F90
	make -f ./makefile all

ifeq ($(FFTLIB),fftw) 
   INCPATHS += -I$(FFTW_HOME)/include
 ifeq ($(PRECISION),double)
   LIBS += -L$(FFTW_HOME)/lib -lfftw3
 else
   LIBS += -L$(FFTW_HOME)/lib -lfftw3f
 endif
endif

FFLAGS += -g 
CFLAGS += -g

ifeq ($(USE_PERFLIB),scalasca)
	FC := scalasca -instrument -user $(FC)

endif

#FORTAN COMPILER FLAGS
FFLAGS += $(MODDIR_FFLAGS) $(DEBUG_WARN)
LIBPATH =

ifeq ($(DEBUG),yes)
	FFLAGS += $(DEBUG_FFLAGS)
	CFLAGS += $(DEBUG_CFLAGS)
else
ifeq ($(DEBUG),valgrind)
	FFLAGS += -O0
	CFLAGS += -O0
else
	FFLAGS += $(OPT_FFLAGS)
	CFLAGS += $(OPT_CFLAGS)
endif
endif

INCPATHS += -I$(OBJDIR) -I$(SRCDIR)

ifeq ($(USE_PERFLIB),hpm)
 HPM_DIR = $(IHPCT_BASE)
 INCPATHS += -I$(HPM_DIR)
 LIBS += -L$(HPM_DIR)/lib -L/afs/ipp/@sys/lib -lhpm_r -lperfhpm_r -lpmapi -lm -lxlsmp
endif
ifeq ($(USE_PERFLIB),perf)
 FFLAGS += $(COMP_WITHPERF)
 LIBS += $(PERFLIB) -lstdc++
endif
ifeq ($(USE_PERFLIB),papi)
 FFLAGS += $(COMP_WITHPAPI)
 INCPATHS += -I$(PAPI_HOME)/include
 LIBS += -L$(PAPI_HOME)/lib -lpapi
endif
ifeq ($(USE_PERFLIB),none)
 OBJLIST += $(OBJDIR)/dummyperf.o
endif

ifeq ($(PRECISION),double)
	FFLAGS += $(DOUBLE_PRECISION_FFLAGS)
endif

LIBS += $(MKL_SCALAPACK_LINKLINE)

all:	tga

include dependencies.mk

TGA_OBJLIST = $(OBJDIR)/test_gyroaverage.o $(LINK_OBJ)

tga:	$(TGA_OBJLIST)
	$(MPFC) -o tga $^ $(LIBS) $(FFLAGS)

$(OBJDIR)/%.o:	$(SRCDIR)/%.F90
	$(MPFC)  $(INCPATHS) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o:	%.F90
	$(MPFC)  $(INCPATHS) $(FFLAGS) -c -o $@ $<

%.s:	%.F90
	$(FC) $(INCPATHS) $(FFLAGS) -S -o $@ $<

%.sc:	%.F90
	$(FC) $(INCPATHS) $(FFLAGS) $<

.PHONY: run
.PHONY:clean
clean:	
	-rm *.o 
	-rm $(OBJDIR)/*.mod
	-rm -f tbm
	-rm -f tfftw
	-rm -f tfftwc





