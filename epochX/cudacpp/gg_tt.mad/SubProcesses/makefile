SHELL := /bin/bash

include ../../Source/make_opts
FFLAGS+= -w

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
CUDACPP_MAKEFILE=cudacpp.mk
# NB1 Using ":=" below instead of "=" is much faster (it only runs the subprocess once instead of many times)
# NB2 Use '|&' in CUDACPP_BUILDDIR to avoid confusing errors about googletest #507
# NB3 Do not add a comment inlined "CUDACPP_BUILDDIR=$(shell ...) # comment" as otherwise a trailing space is included...
# NB4 The variables relevant to the cudacpp Makefile must be explicitly passed to $(shell...)
CUDACPP_MAKEENV:=$(shell echo '$(.VARIABLES)' | tr " " "\n" | egrep "(USEBUILDDIR|AVX|FPTYPE|HELINL|HRDCOD)")
###$(info CUDACPP_MAKEENV=$(CUDACPP_MAKEENV))
###$(info $(foreach v,$(CUDACPP_MAKEENV),$(v)="$($(v))"))
CUDACPP_BUILDDIR:=$(shell $(MAKE) $(foreach v,$(CUDACPP_MAKEENV),$(v)="$($(v))") -f $(CUDACPP_MAKEFILE) -pn |& awk '/Building/{print $$3}' | sed s/BUILDDIR=//)
ifeq ($(CUDACPP_BUILDDIR),)
$(error CUDACPP_BUILDDIR='$(CUDACPP_BUILDDIR)' should not be empty!)
else
$(info CUDACPP_BUILDDIR='$(CUDACPP_BUILDDIR)')
endif
CUDACPP_COMMONLIB=mg5amc_common
CUDACPP_CXXLIB=mg5amc_$(processid_short)_cpp
CUDACPP_CULIB=mg5amc_$(processid_short)_cuda

LIBS = $(LIBDIR)libbias.$(libext) $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libdsample.$(libext) $(LIBDIR)libgeneric.$(libext) $(LIBDIR)libpdf.$(libext) $(LIBDIR)libgammaUPC.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libcernlib.$(libext) $(MADLOOP_LIB) $(LOOP_LIBS)

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
DSIG_cudacpp=driver_cudacpp.o $(patsubst %.f, %_cudacpp.o, $(filter-out auto_dsig.f, $(wildcard auto_dsig*.f)))

SYMMETRY = symmetry.o idenparts.o 

# Binaries

LDFLAGS+=-Wl,--no-relax # avoid 'failed to convert GOTPCREL relocation' error #458

ifeq (,$(wildcard fbridge.inc))
all: $(PROG)_fortran
else
all: $(PROG)_fortran $(CUDACPP_BUILDDIR)/$(PROG)_cpp # also builds $(PROG)_cuda if $(CUDACPP_CULIB) exists (#503)
endif

ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
override OMPFLAGS = -fopenmp
LINKLIBS += -liomp5 # see #578
LINKLIBS += -lintlc # undefined reference to `_intel_fast_memcpy'
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
override OMPFLAGS = -fopenmp
$(CUDACPP_BUILDDIR)/$(PROG)_cpp: LINKLIBS += -L $(shell dirname $(shell $(CXX) -print-file-name=libc++.so)) -lomp # see #604
###else ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
###override OMPFLAGS = -fopenmp # OMP is not supported yet by cudacpp for Apple clang
else
override OMPFLAGS = -fopenmp
endif

$(PROG)_fortran: $(PROCESS) $(DSIG) auto_dsig.o $(LIBS) $(MATRIX) counters.o ompnumthreads.o
	$(FC) -o $(PROG)_fortran $(PROCESS) $(DSIG) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o $(LDFLAGS)

$(LIBS): .libs

.libs: ../../Cards/param_card.dat ../../Cards/run_card.dat
	cd ../../Source; make
	touch $@

ifneq (,$(wildcard fbridge.inc))
$(CUDACPP_BUILDDIR)/.cudacpplibs:
	$(MAKE) -f $(CUDACPP_MAKEFILE)
	touch $@
endif

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override LIBFLAGSRPATH =
else ifeq ($(USEBUILDDIR),1)
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/../$(LIBDIR)/$(CUDACPP_BUILDDIR)'
else
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/$(LIBDIR)'
endif

.PHONY: madevent_fortran_link madevent_cuda_link madevent_cpp_link

madevent_fortran_link: $(PROG)_fortran
	rm -f $(PROG)
	ln -s $(PROG)_fortran $(PROG)

madevent_cpp_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cuda_link: $(CUDACPP_BUILDDIR)/$(PROG)_cuda
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(PROG)

# Building $(PROG)_cpp also builds $(PROG)_cuda if $(CUDACPP_CULIB) exists (improved patch for cpp-only builds #503)
$(CUDACPP_BUILDDIR)/$(PROG)_cpp: $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(LIBS) $(MATRIX) counters.o ompnumthreads.o $(CUDACPP_BUILDDIR)/.cudacpplibs
	$(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) -l$(CUDACPP_CXXLIB) $(LIBFLAGSRPATH) $(LDFLAGS)
	if [ -f $(LIBDIR)/$(CUDACPP_BUILDDIR)/lib$(CUDACPP_CULIB).* ]; then $(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) -l$(CUDACPP_CULIB) $(LIBFLAGSRPATH) $(LDFLAGS); fi

$(CUDACPP_BUILDDIR)/$(PROG)_cuda: $(CUDACPP_BUILDDIR)/$(PROG)_cpp

counters.o: counters.cc timer.h
	$(CXX) -std=c++11 -Wall -Wshadow -Wextra -c $< -o $@

ompnumthreads.o: ompnumthreads.cc ompnumthreads.h
	$(CXX) -I. -std=c++11 -Wall -Wshadow -Wextra $(OMPFLAGS) -c $< -o $@

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(OMPFLAGS)

gensym: $(SYMMETRY) configs.inc $(LIBS)
	$(FC) -o gensym $(SYMMETRY) -L$(LIBDIR) $(LINKLIBS) $(LDFLAGS)

ifeq (,$(wildcard fbridge.inc))
$(LIBDIR)libmodel.$(libext): ../../Cards/param_card.dat
	cd ../../Source/MODEL; make

$(LIBDIR)libgeneric.$(libext): ../../Cards/run_card.dat
	cd ../../Source; make

$(LIBDIR)libpdf.$(libext): 
	cd ../../Source/PDF; make

$(LIBDIR)libgammaUPC.$(libext):
	cd ../../Source/PDF/gammaUPC; make
endif

# Add source so that the compiler finds the DiscreteSampler module.
$(MATRIX): %.o: %.f
	$(FC) $(FFLAGS) $(MATRIX_FLAG) -c $< -I../../Source/ -I../../Source/PDF/gammaUPC
%.o: %.f
	$(FC) $(FFLAGS) -c $< -I../../Source/ -I../../Source/PDF/gammaUPC
%_cudacpp.o: %.f
	$(FC) $(FFLAGS) -c -DMG5AMC_MEEXPORTER_CUDACPP $< -I../../Source/ $(OMPFLAGS) -o $@

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

auto_dsig.o: .libs
driver.o: .libs
driver_cudacpp.o: .libs
$(MATRIX): .libs
genps.o: .libs

# Cudacpp avxall targets

ifneq (,$(wildcard fbridge.inc))

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
avxall: avxnone avxsse4
else ifeq ($(UNAME_P),arm)
avxall: avxnone avxsse4
else
avxall: avxnone avxsse4 avxavx2 avx512y avx512z
endif

avxnone: $(PROG)_fortran $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=none

avxsse4: $(PROG)_fortran $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=sse4

avxavx2: $(PROG)_fortran $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=avx2

avx512y: $(PROG)_fortran $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512y

avx512z: $(PROG)_fortran $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512z

endif

# Clean

clean:
ifeq (,$(wildcard fbridge.inc))
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel
else
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(CUDACPP_BUILDDIR)/$(PROG)_cuda
endif

cleanall:
	make clean
	make -C ../../Source clean
	rm -rf $(LIBDIR)libbias.$(libext)
ifneq (,$(wildcard fbridge.inc))
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleanall
	rm -f $(CUDACPP_BUILDDIR)/.cudacpplibs
endif
	rm -f ../../Source/*.mod ../../Source/*/*.mod
	rm -f .libs

ifneq (,$(wildcard fbridge.inc))
# Target: clean the builds as well as the googletest installation
distclean: cleanall
	$(MAKE) -f $(CUDACPP_MAKEFILE) distclean
endif
