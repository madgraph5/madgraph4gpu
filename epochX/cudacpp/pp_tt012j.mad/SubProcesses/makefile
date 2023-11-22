SHELL := /bin/bash

# Include general setup
OPTIONS_MAKEFILE := ../../Source/make_opts
include $(OPTIONS_MAKEFILE)

# Enable the C preprocessor https://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
MG_FCFLAGS  += -cpp
MG_CXXFLAGS += -I.

all: help cppnative

# Target if user does not specify target
help:
	$(info No target specified.)
	$(info Viable targets are 'cppnative' (default), 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z' and 'cuda')
	$(info Or 'cppall' for all C++ targets)
	$(info Or 'ALL' for all C++ and cuda targets)


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
LIBS = $(LIBDIR)libbias.$(libext) $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libdsample.$(libext) $(LIBDIR)libgeneric.$(libext) $(LIBDIR)libpdf.$(libext) $(LIBDIR)libgammaUPC.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libcernlib.$(libext) $(MADLOOP_LIB) $(LOOP_LIBS)

ifneq ("$(wildcard ../../Source/RUNNING)","")
    LINKLIBS += -lrunning
    LIBS += $(LIBDIR)librunning.$(libext)
endif

SourceDir: ../../Source/*.f ../../Cards/param_card.dat ../../Cards/run_card.dat
	$(MAKE) -C ../../Source

$(LIBS): SourceDir

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
DSIG_cudacpp=driver_cudacpp.o $(patsubst %.f, %_cudacpp.o, $(filter-out auto_dsig.f, $(wildcard auto_dsig*.f)))

SYMMETRY = symmetry.o idenparts.o

# cudacpp targets:
CUDACPP_MAKEFILE := cudacpp.mk
ifneq (,$(wildcard $(CUDACPP_MAKEFILE)))
include $(CUDACPP_MAKEFILE)
endif

ifeq ($(CUDACPP_BUILDDIR),)
$(error CUDACPP_BUILDDIR='$(CUDACPP_BUILDDIR)' should not be empty!)
endif
CUDACPP_COMMONLIB=mg5amc_common
CUDACPP_CXXLIB := $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so
CUDACPP_CULIB  := $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so


# Binaries

$(warning Find a place for the following)
ifeq ($(UNAME),Darwin)
MG_LDFLAGS += -lc++ # avoid 'Undefined symbols' for chrono::steady_clock on macOS (checked with otool -L libmg5amc_gg_ttx_cpp.so)
MG_LDFLAGS += -mmacosx-version-min=11.3 # avoid "ld: warning: object file was built for newer macOS version than being linked"
else
MG_LDFLAGS += -Xlinker --no-relax # avoid 'failed to convert GOTPCREL relocation' error #458 (not supported on macOS)
endif

# Set up OpenMP if supported
OMPFLAGS ?= -fopenmp
ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
LINKLIBS += -liomp5 # see #578
LINKLIBS += -lintlc # undefined reference to `_intel_fast_memcpy'
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
$(CUDACPP_BUILDDIR)/$(PROG)_cpp: LINKLIBS += -L $(shell dirname $(shell $(CXX) -print-file-name=libc++.so)) -lomp # see #604
else ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
OMPFLAGS = # OMP is not supported yet by cudacpp for Apple clang
endif

$(PROG)_fortran: $(PROCESS) $(DSIG) auto_dsig.o $(LIBS) $(MATRIX) counters.o ompnumthreads.o
	$(FC) $(MG_FCFLAGS) $(FCFLAGS) -o $(PROG)_fortran $(PROCESS) $(DSIG) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o $(LDFLAGS)

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifneq ($(UNAME_S),Darwin)
  LIBFLAGSRPATH := -Wl,-rpath,'$$ORIGIN:$$ORIGIN/../$(LIBDIR)/$(CUDACPP_BUILDDIR)'
endif

.PHONY: SourceDir madevent_fortran_link madevent_cuda_link madevent_cpp_link madevent_cppnone_link madevent_cppsse4_link madevent_cppavx2_link madevent_cpp512y_link madevent_cpp512z_link

madevent_fortran_link: $(PROG)_fortran
	rm -f $(PROG)
	ln -s $(PROG)_fortran $(PROG)

madevent_cppnone_link: AVX=none
madevent_cppnone_link: cppnone
	ln -sf $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cppavx2_link: AVX=avx2
madevent_cppavx2_link: cppavx2
	ln -sf $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cpp512y_link: AVX=512y
madevent_cpp512y_link: cppavx512y
	ln -sf $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cpp512z_link: AVX=512z
madevent_cpp512z_link: cppavx512z
	ln -sf $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cuda_link: AVX=cuda
madevent_cuda_link: cuda
	ln -sf $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(PROG)

$(CUDACPP_BUILDDIR)/$(PROG)_cpp: $(LIBS) $(CUDACPP_CXXLIB) $(PROCESS) $(DSIG_cudacpp) auto_dsig.o  $(MATRIX) counters.o ompnumthreads.o
	$(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cpp -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CXXLIB) $(LIBFLAGSRPATH) $(MG_LDFLAGS) $(LDFLAGS)

$(CUDACPP_BUILDDIR)/$(PROG)_cuda: $(LIBS) $(CUDACPP_CULIB) $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) counters.o ompnumthreads.o
	$(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cuda -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CULIB) $(LIBFLAGSRPATH) $(MG_LDFLAGS) $(LDFLAGS)

counters.o: counters.cc timer.h
	$(CXX) $(MG_CXXFLAGS) $(CXXFLAGS) -c $< -o $@

ompnumthreads.o: ompnumthreads.cc ompnumthreads.h
	$(CXX) $(MG_CXXFLAGS) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(OMPFLAGS)

gensym: $(SYMMETRY) configs.inc $(LIBS)
	$(FC) -o gensym $(SYMMETRY) -L$(LIBDIR) $(LINKLIBS) $(LDFLAGS)


# Add source so that the compiler finds the DiscreteSampler module.
$(MATRIX): %.o: %.f
	$(FC) $(MG_FCFLAGS) $(FCFLAGS) $(MATRIX_FLAG) -c $< -I../../Source/ -I../../Source/PDF/gammaUPC
%.o $(CUDACPP_BUILDDIR)/%.o: %.f
	$(FC) $(MG_FCFLAGS) $(FCFLAGS) -I../../Source/ -I../../Source/PDF/gammaUPC -c $< -o $@
%_cudacpp.o: %.f
	$(FC) $(MG_FCFLAGS) $(FCFLAGS) -c -DMG5AMC_MEEXPORTER_CUDACPP $< -I../../Source/ $(OMPFLAGS) -o $@

# Dependencies

driver.f: genps.inc
symmetry.o: genps.inc nexternal.inc configs.inc run_config.inc ../../Source/run_card.inc
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
../../Source/discretesampler.mod: ../../Source/DiscreteSampler.f
	$(MAKE) -C ../../Source

auto_dsig.o: ../../Source/discretesampler.mod
driver.o: ../../Source/discretesampler.mod
driver_cudacpp.o: ../../Source/discretesampler.mod
$(MATRIX): ../../Source/discretesampler.mod
genps.o: ../../Source/discretesampler.mod

# Cudacpp avxall targets

ifeq ($(UNAME_P),ppc64le)
cppall: cppnative cppnone cppsse4
else ifeq ($(UNAME_P),arm)
cppall: cppnative cppnone cppsse4
else
cppall: cppnative cppnone cppsse4 cppavx2 cppavx512y cppavx512z
endif

ALL: cppall cuda

# Clean all architecture-specific builds:
clean:
	$(RM) *.o gensym $(PROG) $(PROG)_*
	$(RM) -rf build.*/*{.o,.so,.exe,.dylib,madevent_*}

# Clean one architecture-specific build
clean%:
	$(RM) -r build.$*_*/*{.o,.so,.exe,.dylib,madevent_*}

# Clean common source directories (interferes with other P*)
cleansrc:
	make -C ../../Source clean
	$(RM) ../../Source/*.mod ../../Source/*/*.mod
	$(RM) -r $(LIBDIR)libbias.$(libext)
	if [ -d ../../src ]; then $(MAKE) -C ../../src -f cudacpp_src.mk clean; fi
