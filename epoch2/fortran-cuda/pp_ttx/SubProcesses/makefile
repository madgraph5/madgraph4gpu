include ../../Source/make_opts
FFLAGS+= -w

EPOCHDIR=epochX
LANGDIR=cudacpp
PROCESSDIR=gg_tt.auto
MEBASEDIR=../../../../../$(EPOCHDIR)/$(LANGDIR)/$(PROCESSDIR)/
MELIBDIR=$(MEBASEDIR)/lib
MEINCDIR=$(MEBASEDIR)/include

CUINCDIR=/usr/local/cuda-11.4/include/
TOOLSINCDIR=../../../../../tools

GCCHOME=/opt/rh/gcc-toolset-10/root/usr
GCCBIN=$(GCCHOME)/bin
GPPEXE=$(GCCBIN)/g++

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
    CUDA_HOME  = $(patsubst %/bin/nvcc,%,$(NVCC))
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

# Set the default AVX (vectorization) choice
ifeq ($(AVX),)
  ifeq ($(UNAME_P),ppc64le)
    ###override AVX = none
    override AVX = sse4
  else ifeq ($(shell grep -m1 -c avx512vl /proc/cpuinfo)$(shell $(CXX) --version | grep ^clang),1)
    override AVX = 512y
    ###$(info Using AVX='$(AVX)' as no user input exists)
  else
    override AVX = avx2
    ifneq ($(shell grep -m1 -c avx512vl /proc/cpuinfo),1)
      $(warning Using AVX='$(AVX)' because host does not support avx512vl)
    else
      $(warning Using AVX='$(AVX)' because this is faster than avx512vl for clang)
    endif
  endif
else
  ###$(info Using AVX='$(AVX)' according to user input)
endif

# Set the build flags appropriate to each AVX choice (example: "make AVX=none")
# [NB MGONGPU_PVW512 is needed because "-mprefer-vector-width=256" is not exposed in a macro]
# [See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96476]
$(info AVX=$(AVX))
ifeq ($(UNAME_P),ppc64le)
  ifeq ($(AVX),sse4)
    override AVXFLAGS = -D__SSE4_2__ # Power9 VSX with 128 width (VSR registers)
  else ifneq ($(AVX),none)
    $(error Unknown AVX='$(AVX)': only 'none' and 'sse4' are supported on PowerPC for the moment)
  endif
else
  ifeq ($(AVX),sse4)
    override AVXFLAGS = -march=nehalem # SSE4.2 with 128 width (xmm registers)
  else ifeq ($(AVX),avx2)
    override AVXFLAGS = -march=haswell # AVX2 with 256 width (ymm registers) [DEFAULT for clang]
  else ifeq ($(AVX),512y)
    override AVXFLAGS = -march=skylake-avx512 -mprefer-vector-width=256 # AVX512 with 256 width (ymm registers) [DEFAULT for gcc]
  else ifeq ($(AVX),512z)
    override AVXFLAGS = -march=skylake-avx512 -DMGONGPU_PVW512 # AVX512 with 512 width (zmm registers)
  else ifneq ($(AVX),none)
    $(error Unknown AVX='$(AVX)': only 'none', 'sse4', 'avx2', '512y' and '512z' are supported)
  endif
endif
# For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations?
CXXFLAGS+= $(AVXFLAGS)




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

BRIDGE = f2cu_bridge.o f2cpp_bridge.o

MEOBJS=$(MELIBDIR)/gCPPProcess.o $(MELIBDIR)/CPPProcess.o $(MELIBDIR)/libmodel_sm.a

CUDA =

# Binaries

# roiser@itscrd02:P1_gg_ttx (fortran-cpp-workflow)$ /usr/local/cuda-11.4/bin/nvcc -o madevent  driver.o myamp.o genps.o unwgt.o setcuts.o get_color.o cuts.o cluster.o reweight.o initcluster.o addmothers.o setscales.o idenparts.o dummy_fct.o auto_dsig1.o auto_dsig.o  matrix1_optim.o fbridge.o -Xcompiler -fopenmp -L../../lib/ -L /usr/local/cuda-11.4/lib64/ -ldhelas -ldsample -lmodel -lgeneric -lpdf -lcernlib  -lbias -lstdc++ -lcudart ./golden_epochX4/gCPPProcess.o ../../Source/lib/libmodel_sm.a -lgfortran



$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) $(BRIDGE) $(CUDA) $(MEOBJS)
	$(NVCC) -o $(PROG) $(PROCESS) $(MATRIX) $(BRIDGE) $(CUDA) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(CULIBS) $(MEOBJS) -Xcompiler -fopenmp -lgfortran -lcudart

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL) $(BRIDGE) $(CUDA) $(MEOBJS)
	$(nvcc) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(BRIDGE) $(CUDA) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(CULIBS) $(MEOBJS) -Xcompiler -fopenmp -lgfortran -lcudart

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
%.o: %.cpp bridge.h
	$(CXX) $(CXXFLAGS) -O3 -g -std=c++17 -lineinfo -I. -I$(MEINCDIR) -I$(TOOLSINCDIR) -c $< -o $@
%.o: %.cu bridge.h
	$(NVCC) $(CUFLAGS) -O3 -g -std=c++17 -lineinfo -use_fast_math -I. -I$(MEINCDIR) -I$(CUINCDIR) -I$(TOOLSINCDIR) -ccbin $(GPPEXE) -arch=compute_70 -DUSE_NVTX -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_CURAND_ONDEVICE -c $< -o $@


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
	cat G1/input2_app.txt | ./$(PROG)
