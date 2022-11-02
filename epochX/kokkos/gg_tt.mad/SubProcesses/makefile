include ../../Source/make_opts
FFLAGS+= -w

# Add Kokkos flags and variables
KOKKOSPATH_CUDA ?= $(KOKKOS_HOME)
KOKKOSPATH_OMP ?= $(KOKKOS_HOME)
KOKKOSPATH_INTEL ?= $(KOKKOS_HOME)
KOKKOSPATH_HIP ?= $(KOKKOS_HOME)

INCDIR=../../src
CXXFLAGS=-O3 -ffast-math -I$(INCDIR) --std=c++17

LDFLAGS=-lstdc++fs

CUDA_ARCH_NUM ?= 70
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_CUDA)/lib64/libkokkoscore.a -lnvToolsExt --openmp

CXX ?= g++
ICPX ?= icpx
CLANG ?= clang++
HIPCC ?= hipcc

OPENMP_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_OMP)/lib64/libkokkoscore.a -ldl --openmp

INTEL_BACKEND ?= gen9
INTEL_SYCL ?= 0
ifeq ($(INTEL_SYCL),1)
# build for SYCL
INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a $(KOKKOSPATH_INTEL)/lib64/libkokkoscontainers.a -fiopenmp -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend "-device ${INTEL_BACKEND}"
else
# build for OPENMP
INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a $(KOKKOSPATH_INTEL)/lib64/libkokkoscontainers.a -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device ${INTEL_BACKEND}"
endif

HIP_ARCH_NUM ?= gfx906
HIP_LDFLAGS=$(LDFLAGS) -L $(KOKKOSPATH_HIP)/lib64  -lkokkoscore -ldl -fopenmp 

# Enable the C preprocessor https://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
FFLAGS+= -cpp 

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif
ifeq ($(USECCACHE)$(shell echo $(FC) | grep ccache),1)
  override FC:=ccache $(FC)
endif

# Load additional dependencies of the bias module, if present
ifeq (,$(wildcard ../bias_dependencies))
BIASDEPENDENCIES =
else
include ../bias_dependencies
endif

# Definitions

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

LINKLIBS = $(LINK_MADLOOP_LIB) $(LINK_LOOP_LIBS) -L$(LIBDIR) -ldhelas -ldsample -lmodel -lgeneric -lpdf -lcernlib $(llhapdf) -lbias 

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')
KOKKOS_MAKEFILE=kokkos.mk
KOKKOSSRC_MAKEFILE=kokkos_src.mk
# NB1 Using ":=" below instead of "=" is much faster (it only runs the subprocess once instead of many times)
# NB2 Do not add a comment inlined "KOKKOS_BUILDDIR=$(shell ...) # comment" as otherwise a trailing space is included...
# NB3 The variables relevant to the kokkos Makefile must be explicitly passed to $(shell...)
KOKKOS_MAKEENV:=$(shell echo '$(.VARIABLES)' | tr " " "\n" | egrep "(USEBUILDDIR|AVX|FPTYPE|HELINL|HRDCOD)")
###$(info KOKKOS_MAKEENV=$(KOKKOS_MAKEENV))
###$(info $(foreach v,$(KOKKOS_MAKEENV),$(v)="$($(v))"))
#KOKKOS_BUILDDIR:=$(shell $(MAKE) $(foreach v,$(KOKKOS_MAKEENV),$(v)="$($(v))") -f $(KOKKOS_MAKEFILE) -pn | awk '/Building/{print $$3}' | sed s/BUILDDIR=//)
KOKKOS_BUILDDIR:=build
###$(info KOKKOS_BUILDDIR='$(KOKKOS_BUILDDIR)')
KOKKOS_COMMONLIB=mg5amc_common
KOKKOS_BRIDGELIB=mg5amc_$(processid_short)_kokkos

LIBS = $(LIBDIR)libbias.$(libext) $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libdsample.$(libext) $(LIBDIR)libgeneric.$(libext) $(LIBDIR)libpdf.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libcernlib.$(libext) $(MADLOOP_LIB) $(LOOP_LIBS)

ifneq ("$(wildcard ../../Source/RUNNING)","")
    LINKLIBS += -lrunning
    LIBS += $(LIBDIR)librunning.$(libext) 
endif


# Source files

MATRIX_HEL = $(patsubst %.f,%.o,$(wildcard matrix*_orig.f))
MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*_optim.f))
ifeq ($(strip $(MATRIX_HEL)),)
        MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*.f))
endif


PROCESS= myamp.o genps.o unwgt.o setcuts.o get_color.o \
         cuts.o cluster.o reweight.o initcluster.o addmothers.o setscales.o \
	 idenparts.o dummy_fct.o

DSIG=driver.o $(patsubst %.f, %.o, $(filter-out auto_dsig.f, $(wildcard auto_dsig*.f)))
DSIG_kokkos=driver_kokkos.o $(patsubst %.f, %_kokkos.o, $(filter-out auto_dsig.f, $(wildcard auto_dsig*.f)))

SYMMETRY = symmetry.o idenparts.o 

# Binaries

LDFLAGS+=-Wl,--no-relax # avoid 'failed to convert GOTPCREL relocation' error #458

ifeq (,$(wildcard fbridge.inc))
all: $(PROG)
else
all: $(PROG)
endif

$(PROG): $(PROCESS) $(DSIG) auto_dsig.o $(LIBS) $(MATRIX) counters.o
	$(FC) -o $(PROG) $(PROCESS) $(DSIG) auto_dsig.o $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o -lstdc++

$(LIBS): .libs

.libs: ../../Cards/param_card.dat ../../Cards/run_card.dat
	cd ../../Source; make
	touch $@

ifneq (,$(wildcard fbridge.inc))
$(KOKKOS_BUILDDIR)/.kokkoslibs_cuda:
	$(MAKE) -f $(KOKKOS_MAKEFILE) ../../lib/libmodel_sm_bridge_cuda.a
	$(MAKE) -C ../../src -f $(KOKKOSSRC_MAKEFILE) ../lib/libmodel_sm_cuda.a
	touch $@

$(KOKKOS_BUILDDIR)/.kokkoslibs_openmp:
	$(MAKE) -f $(KOKKOS_MAKEFILE) ../../lib/libmodel_sm_bridge_openmp.a
	$(MAKE) -C ../../src -f $(KOKKOSSRC_MAKEFILE) ../lib/libmodel_sm_openmp.a
	touch $@

$(KOKKOS_BUILDDIR)/.kokkoslibs_intel:
	$(MAKE) -f $(KOKKOS_MAKEFILE) ../../lib/libmodel_sm_bridge_intel.a
	$(MAKE) -C ../../src -f $(KOKKOSSRC_MAKEFILE) ../lib/libmodel_sm_intel.a
	touch $@

$(KOKKOS_BUILDDIR)/.kokkoslibs_hip:
	$(MAKE) -f $(KOKKOS_MAKEFILE) ../../lib/libmodel_sm_bridge_hip.a
	$(MAKE) -C ../../src -f $(KOKKOSSRC_MAKEFILE) ../lib/libmodel_sm_hip.a
	touch $@
endif

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override LIBFLAGSRPATH =
else ifeq ($(USEBUILDDIR),1)
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/../$(LIBDIR)/$(KOKKOS_BUILDDIR)'
else
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/$(LIBDIR)'
endif

cuda: $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(KOKKOS_BUILDDIR)/.kokkoslibs_cuda
	$(CXX) $(KOKKOSFLAGS) -o $(KOKKOS_BUILDDIR)/g$(PROG)_kokkos $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(MATRIX) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a $(LIBDIR)/lib$(KOKKOS_COMMONLIB).a $(LINKLIBS) $(LDFLAGS) $(CUDA_LDFLAGS) $(BIASDEPENDENCIES) counters.o -lgomp -lgfortran $(LIBFLAGSRPATH) -lstdc++fs

openmp: $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(KOKKOS_BUILDDIR)/.kokkoslibs_openmp
	$(CXX) $(KOKKOSFLAGS) -o $(KOKKOS_BUILDDIR)/g$(PROG)_kokkos $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(MATRIX) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a $(LIBDIR)/lib$(KOKKOS_COMMONLIB).a $(LINKLIBS) $(LDFLAGS) $(OPENMP_LDFLAGS) $(BIASDEPENDENCIES) counters.o -lgomp -lgfortran $(LIBFLAGSRPATH) -lstdc++fs

intel: $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(KOKKOS_BUILDDIR)/.kokkoslibs_intel
	$(CXX) $(KOKKOSFLAGS) -o $(KOKKOS_BUILDDIR)/g$(PROG)_kokkos $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(MATRIX) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a $(LIBDIR)/lib$(KOKKOS_COMMONLIB).a $(LINKLIBS) $(LDFLAGS) $(INTEL_LDFLAGS) $(BIASDEPENDENCIES) counters.o -lgomp -lgfortran $(LIBFLAGSRPATH) -lstdc++fs

hip: $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(KOKKOS_BUILDDIR)/.kokkoslibs_hip
	$(CXX) $(KOKKOSFLAGS) -o $(KOKKOS_BUILDDIR)/g$(PROG)_kokkos $(PROCESS) $(DSIG_kokkos) auto_dsig.o $(MATRIX) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a $(LIBDIR)/lib$(KOKKOS_COMMONLIB).a $(LINKLIBS) $(LDFLAGS) $(HIP_LDFLAGS) $(BIASDEPENDENCIES) counters.o -lgomp -lgfortran $(LIBFLAGSRPATH) -lstdc++fs

counters.o: counters.cpp timer.h
	$(CXX) -std=c++11 -Wall -Wshadow -Wextra -c $< -o $@

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp

gensym: $(SYMMETRY) configs.inc $(LIBS)
	$(FC) -o gensym $(SYMMETRY) -L$(LIBDIR) $(LINKLIBS) $(LDFLAGS)

ifeq (,$(wildcard fbridge.inc))
$(LIBDIR)libmodel.$(libext): ../../Cards/param_card.dat
	cd ../../Source/MODEL; make

$(LIBDIR)libgeneric.$(libext): ../../Cards/run_card.dat
	cd ../../Source; make

$(LIBDIR)libpdf.$(libext): 
	cd ../../Source/PDF; make
endif

# Add source so that the compiler finds the DiscreteSampler module.
$(MATRIX): %.o: %.f
	$(FC) $(FFLAGS) $(MATRIX_FLAG) -c $< -I../../Source/ -fopenmp
%.o: %.f
	$(FC) $(FFLAGS) -c $< -I../../Source/ -fopenmp -o $@
%_kokkos.o: %.f
	$(FC) $(FFLAGS) -c -DMG5AMC_MEEXPORTER_KOKKOS $< -I../../Source/ -fopenmp -o $@

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

# Extra dependencies on discretesampler.mod

auto_dsig.o: .libs
driver.o: .libs
driver_kokkos.o: .libs
$(MATRIX): .libs
genps.o: .libs

# kokkos avxall targets

ifneq (,$(wildcard fbridge.inc))

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
avxall: avxnone avxsse4
else ifeq ($(UNAME_P),arm)
avxall: avxnone avxsse4
else
avxall: avxnone avxsse4 avxavx2 avx512y avx512z
endif

avxnone: $(PROG) $(DSIG_kokkos)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=none

avxsse4: $(PROG) $(DSIG_kokkos)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=sse4

avxavx2: $(PROG) $(DSIG_kokkos)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=avx2

avx512y: $(PROG) $(DSIG_kokkos)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512y

avx512z: $(PROG) $(DSIG_kokkos)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512z

endif

# Clean

clean:
ifeq (,$(wildcard fbridge.inc))
	$(RM) *.o gensym $(PROG) $(PROG)_forhel
else
	$(RM) *.o gensym $(PROG) $(PROG)_forhel $(KOKKOS_BUILDDIR)/*$(PROG)_kokkos
endif

cleanall:
	make clean
	make -C ../../Source clean
	rm -rf $(LIBDIR)libbias.$(libext)
ifneq (,$(wildcard fbridge.inc))
	$(MAKE) -f $(KOKKOS_MAKEFILE) cleanall
	rm -f $(KOKKOS_BUILDDIR)/.kokkoslibs
endif
	rm -f ../../Source/*.mod ../../Source/*/*.mod
	rm -f .libs
