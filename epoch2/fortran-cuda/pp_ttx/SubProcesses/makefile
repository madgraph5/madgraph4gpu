include ../../Source/make_opts
FFLAGS+= -w

# Load additional dependencies of the bias module, if present
ifeq (,$(wildcard ../bias_dependencies))
BIASDEPENDENCIES =
else
include ../bias_dependencies
endif

# Definitions
ifndef CUDA_HOME
  NVCC ?= $(shell which nvcc 2>/dev/null)
  ifneq ($(NVCC),)
    CUDA_HOME  = $(patsubst %bin/nvcc,%,$(NVCC))
    CUDA_HOME := $(warning No CUDA_HOME exported. Using "$(CUDA_HOME)") $(CUDA_HOME)
  endif
endif

NVCC = $(CUDA_HOME)/bin/nvcc
LIBDIR = ../../lib/
BINDIR = ../../bin/
PROG   = madevent

ifneq ("$(wildcard ../MadLoop_makefile_definitions)","")
    include ../MadLoop_makefile_definitions
else
    LINK_LOOP_LIBS =
    LOOP_LIBS =
    LOOP_INCLUDE =
    LINK_MADLOOP_LIB =
    MADLOOP_LIB =
endif

LINKLIBS = $(LINK_MADLOOP_LIB) $(LINK_LOOP_LIBS) -L../../lib/ -ldhelas -ldsample -lmodel -lgeneric -lpdf -lcernlib $(llhapdf) -lbias

LIBS = $(LIBDIR)libbias.$(libext) $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libdsample.$(libext) $(LIBDIR)libgeneric.$(libext) $(LIBDIR)libpdf.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libcernlib.$(libext) $(MADLOOP_LIB) $(LOOP_LIBS)
CULIBS = -L$(CUDA_HOME)/lib64/ -lcudart

# Source files

MATRIX_HEL = $(patsubst %.f,%.o,$(wildcard matrix*_orig.f))
MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*_optim.f))
ifeq ($(strip $(MATRIX_HEL)),)
        MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*.f))
endif


PROCESS= driver.o myamp.o genps.o unwgt.o setcuts.o get_color.o \
         cuts.o cluster.o reweight.o initcluster.o addmothers.o setscales.o \
	 idenparts.o dummy_fct.o \
         $(patsubst %.f,%.o,$(wildcard auto_dsig*.f)) \

SYMMETRY = symmetry.o idenparts.o

BRIDGE = fbridge.o

CUDA = transpose.o

# Binaries

$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) $(BRIDGE)
	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(BRIDGE) $(CUDA) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(CULIBS) -fopenmp

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL) $(BRIDGE)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(BRIDGE) $(CUDA) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(CULIBS) -fopenmp

gensym: $(SYMMETRY) configs.inc $(LIBDIR)libmodel.$(libext) $(LIBDIR)libgeneric.$(libext)
	$(FC) -o gensym $(SYMMETRY) -L../../lib/ -lmodel -lgeneric -lpdf $(llhapdf) $(LDFLAGS)

$(LIBDIR)libmodel.$(libext): ../../Cards/param_card.dat
	cd ../../Source/MODEL; make

$(LIBDIR)libgeneric.$(libext): ../../Cards/run_card.dat
	cd ../../Source; make

$(LIBDIR)libpdf.$(libext):
	cd ../../Source/PDF; make

# Add source so that the compiler finds the DiscreteSampler module.
$(MATRIX): %.o: %.f
	$(FC) $(FFLAGS) $(MATRIX_FLAG) -c $< -I../../Source/ -fopenmp
%.o: %.f
	$(FC) $(FFLAGS) -c $< -I../../Source/ -fopenmp
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
%.o: %.cu
	$(NVCC) $(CUFLAGS) -c $< -o $@


# Dependencies

driver.f: genps.inc
symmetry.o: genps.inc nexternal.inc configs.inc run_config.inc
genps.o: genps.inc nexternal.inc configs.inc
dummy_fct.0: run.inc genps.inc
cuts.o: genps.inc nexternal.inc pmass.inc
setcuts.o: genps.inc run_config.inc
invarients.o: genps.inc nexternal.inc
myamp.o: props.inc genps.inc nexternal.inc
reweight.o: sudakov.inc cluster.inc sudakov.inc run.inc message.inc \
	    run_config.inc
cluster.o: cluster.inc genps.inc nexternal.inc message.inc
addmothers.o: genps.inc nexternal.inc symswap.inc message.inc
unwgt.o: genps.inc nexternal.inc symswap.inc cluster.inc run.inc message.inc \
	 run_config.inc
initcluster.o: message.inc

clean:
	$(RM) *.o gensym madevent madevent_forhel

run: run1

run1: $(PROG)
	cat G1/input_app.txt | ./$(PROG)
