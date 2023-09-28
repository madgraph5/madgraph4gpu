SHELL := /bin/bash

include ../../Source/make_opts
FFLAGS+= -w

# Enable the C preprocessor https://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
FFLAGS+= -cpp

# Compile counters with -O3 as in the cudacpp makefile (avoid being "unfair" to Fortran #740)
CXXFLAGS = -O3 -Wall -Wshadow -Wextra

# Sets correct target based on MAKECMDGOALS
ifeq ($(MAKECMDGOALS),)
    TARGET := undefined
else
    ifneq (,$(findstring madevent_*,$(MAKECMDGOALS)))
        TARGET := $(patsubst %_%,%,$(MAKECMDGOALS))
    else
        TARGET := $(MAKECMDGOALS)
    endif
endif

# Default goal
.DEFAULT_GOAL := usage

# Target if user does not specify target
usage:
	$(error Unknown target='$(TARGET)': only 'cppnone', 'cppsse4', 'cppavx2', 'cpp512y', 'cpp512z' and 'cuda' are supported!)

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
CUDACPP_MAKEENV:=$(shell echo '$(.VARIABLES)' | tr " " "\n" | egrep "(USEBUILDDIR|FPTYPE|HELINL|HRDCOD)")
###$(info CUDACPP_MAKEENV=$(CUDACPP_MAKEENV))
###$(info $(foreach v,$(CUDACPP_MAKEENV),$(v)="$($(v))"))
CUDACPP_BUILDDIR:=$(shell $(MAKE) $(foreach v,$(CUDACPP_MAKEENV),$(v)="$($(v))") -f $(CUDACPP_MAKEFILE) -pn $(TARGET) 2>/dev/null |& awk '/Building/{print $$3}' | sed s/BUILDDIR=//)
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

all: $(PROG)_fortran $(CUDACPP_BUILDDIR)/$(PROG)_cpp # also builds $(PROG)_cuda if $(CUDACPP_CULIB) exists (#503)

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

.PHONY: madevent_fortran_link madevent_cuda_link madevent_cpp_link madevent_cppnone_link madevent_cppsse4_link madevent_cppavx2_link madevent_cpp512y_link madevent_cpp512z_link

madevent_fortran_link: $(PROG)_fortran
	rm -f $(PROG)
	ln -s $(PROG)_fortran $(PROG)

# Only used for testing, should be removed in the future
madevent_cpp_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cppavx2
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cppnone_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cppnone
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cppsse4_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cppsse4
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cppavx2_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cppavx2
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cpp512y_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cpp512y
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cpp512z_link: $(CUDACPP_BUILDDIR)/$(PROG)_cpp cpp512z
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROG)

madevent_cuda_link: $(CUDACPP_BUILDDIR)/$(PROG)_cuda cuda
	rm -f $(PROG)
	ln -s $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(PROG)

# Building $(PROG)_cpp also builds $(PROG)_cuda if $(CUDACPP_CULIB) exists (improved patch for cpp-only builds #503)
$(CUDACPP_BUILDDIR)/$(PROG)_cpp: $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(LIBS) $(MATRIX) counters.o ompnumthreads.o
	$(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) -l$(CUDACPP_CXXLIB) $(LIBFLAGSRPATH) $(LDFLAGS)

$(CUDACPP_BUILDDIR)/$(PROG)_cuda: $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(LIBS) $(MATRIX) counters.o ompnumthreads.o
	$(FC) -o $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(MATRIX) $(LINKLIBS) $(BIASDEPENDENCIES) $(OMPFLAGS) counters.o ompnumthreads.o -L$(LIBDIR)/$(CUDACPP_BUILDDIR) -l$(CUDACPP_COMMONLIB) -l$(CUDACPP_CULIB) $(LIBFLAGSRPATH) $(LDFLAGS)

counters.o: counters.cc timer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

ompnumthreads.o: ompnumthreads.cc ompnumthreads.h
	$(CXX) -I. $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) $(OMPFLAGS)

gensym: $(SYMMETRY) configs.inc $(LIBS)
	$(FC) -o gensym $(SYMMETRY) -L$(LIBDIR) $(LINKLIBS) $(LDFLAGS)

###ifeq (,$(wildcard fbridge.inc)) # Pointless: fbridge.inc always exists as this is the cudacpp-modified makefile!
###$(LIBDIR)libmodel.$(libext): ../../Cards/param_card.dat
###	cd ../../Source/MODEL; make
###
###$(LIBDIR)libgeneric.$(libext): ../../Cards/run_card.dat
###	cd ../../Source; make
###
###$(LIBDIR)libpdf.$(libext): 
###	cd ../../Source/PDF; make
###
###$(LIBDIR)libgammaUPC.$(libext):
###	cd ../../Source/PDF/gammaUPC; make
###endif

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

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
cppall:
	$(MAKE) USEBUILDDIR=1 cppnone
	$(MAKE) USEBUILDDIR=1 cppsse4
else ifeq ($(UNAME_P),arm)
cppall:
	$(MAKE) USEBUILDDIR=1 cppnone
	$(MAKE) USEBUILDDIR=1 cppsse4
else
cppall:
	$(MAKE) USEBUILDDIR=1 cppnone
	$(MAKE) USEBUILDDIR=1 cppsse4
	$(MAKE) USEBUILDDIR=1 cppavx2
	$(MAKE) USEBUILDDIR=1 cpp512y
	$(MAKE) USEBUILDDIR=1 cpp512z
endif

ALL:
	$(MAKE) USEBUILDDIR=1 cppnone
	$(MAKE) USEBUILDDIR=1 cppsse4
	$(MAKE) USEBUILDDIR=1 cppavx2
	$(MAKE) USEBUILDDIR=1 cpp512y
	$(MAKE) USEBUILDDIR=1 cpp512z
	$(MAKE) USEBUILDDIR=1 cuda

cppnone:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cppnone

cppsse4:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cppsse4

cppavx2:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cppavx2

cpp512y:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cpp512y

cpp512z:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cpp512z

cuda:
	$(MAKE) -f $(CUDACPP_MAKEFILE) cuda

# Clean (NB: 'make clean' in Source calls 'make clean' in all P*)

BUILD_DIRS := $(wildcard build.*)
BUILD_TARGET_DIR := $(shell echo $(BUILD_DIRS) | awk -F '[._]' '{print $$2}')
NUM_BUILD_DIRS := $(words $(BUILD_DIRS))

clean:
ifeq ($(USEBUILDDIR),1)
ifeq ($(NUM_BUILD_DIRS),1)
	$(info USEBUILDDIR=1, only cleaning build dir.)
	$(RM) -r $(BUILD_DIRS)
	$(MAKE) -f $(CUDACPP_MAKEFILE) clean$(BUILD_TARGET_DIR)
else ifeq ($(NUM_BUILD_DIRS),0)
	$(error USEBUILDDIR=1, but no build directories are found.)
else
	$(error Multiple BUILDDIR's found! Use 'cleannone', 'cleansse4', 'cleanavx2', 'clean512y','clean512z', 'cleancuda' or 'cleanall'.)
endif
else
	$(info USEBUILDDIR=0 only cleaning build with no dir.)
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel ./$(PROG)_*
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleandir
	$(RM) $(CUDACPP_BUILDDIR)/.*libs
endif

cleannone: # Clean builds: fortran in this Pn; cpp executables for one AVX in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.none_*/$(PROG)_cpp
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleannone

cleansse4: # Clean build: fortran in this Pn; cpp executables for cppsse4 in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.sse4_*/$(PROG)_cpp
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleansse4

cleanavx2: # Clean build: fortran in this Pn; cpp executables for cppavx2 in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.avx2_*/$(PROG)_cpp
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleanavx2

clean512y: # Clean build: fortran in this Pn; cpp executables for cpp512y in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.512y_*/$(PROG)_cpp
	$(MAKE) -f $(CUDACPP_MAKEFILE) clean512y

clean512z: # Clean build: fortran in this Pn; cpp executables for cpp512z in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.512z_*/$(PROG)_cpp
	$(MAKE) -f $(CUDACPP_MAKEFILE) clean512z

cleancuda: # Clean build: fortran in this Pn; cuda executables in this Pn
	$(RM) *.o gensym $(PROG) $(PROG)_fortran $(PROG)_forhel build.cuda_*/$(PROG)_cuda
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleancuda

cleanavxs: clean # Clean builds: fortran in this Pn; cudacpp for all AVX in this Pn and in src
	$(MAKE) -f $(CUDACPP_MAKEFILE) cleanall
	$(RM) $(CUDACPP_BUILDDIR)/.*libs
	$(RM) .libs

cleanall: # Clean builds: fortran in all P* and in Source; cudacpp for all AVX in all P* and in src
	make -C ../../Source cleanall
	$(RM) -r $(LIBDIR)libbias.$(libext)
	$(RM) ../../Source/*.mod ../../Source/*/*.mod

distclean: cleanall # Clean all fortran and cudacpp builds as well as the googletest installation
	$(MAKE) -f $(CUDACPP_MAKEFILE) distclean
