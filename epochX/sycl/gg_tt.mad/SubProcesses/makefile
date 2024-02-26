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
SYCL_MAKEFILE=sycl.mk
# NB1 Using ":=" below instead of "=" is much faster (it only runs the subprocess once instead of many times)
# NB2 Use '|&' in SYCL_BUILDDIR to avoid confusing errors about googletest #507
# NB3 Do not add a comment inlined "SYCL_BUILDDIR=$(shell ...) # comment" as otherwise a trailing space is included...
# NB4 The variables relevant to the sycl Makefile must be explicitly passed to $(shell...)
SYCL_MAKEENV:=$(shell echo '$(.VARIABLES)' | tr " " "\n" | egrep "(USEBUILDDIR|FPTYPE|HELINL|HRDCOD|NTPBMAX|CXTYPE|VECLVL)")
###$(info SYCL_MAKEENV=$(SYCL_MAKEENV))
###$(info $(foreach v,$(SYCL_MAKEENV),$(v)="$($(v))"))
SYCL_BUILDDIR:=$(shell $(MAKE) $(foreach v,$(SYCL_MAKEENV),$(v)="$($(v))") -f $(SYCL_MAKEFILE) -pn |& awk '/Building/{print $$3}' | sed s/BUILDDIR=//)
ifeq ($(SYCL_BUILDDIR),)
$(error SYCL_BUILDDIR='$(SYCL_BUILDDIR)' should not be empty!)
else
$(info SYCL_BUILDDIR='$(SYCL_BUILDDIR)')
endif
SYCL_COMMONLIB=mg5amc_common
SYCL_SYCLLIB=mg5amc_$(processid_short)_sycl

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
DSIG_sycl=driver_sycl.o $(patsubst %.f, %_sycl.o, $(filter-out auto_dsig.f, $(wildcard auto_dsig*.f)))

SYMMETRY = symmetry.o idenparts.o 

# Binaries

LDFLAGS+=-Wl,--no-relax # avoid 'failed to convert GOTPCREL relocation' error #458

ifeq (,$(wildcard fbridge.inc))
all: $(PROG)
else
all: $(PROG) $(SYCL_BUILDDIR)/$(PROG)_sycl
endif


$(PROG): $(PROCESS) $(DSIG) auto_dsig.o $(LIBS) $(MATRIX) counters.o
	$(FC) -o $(PROG) $(PROCESS) $(DSIG) auto_dsig.o $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o

$(LIBS): .libs

.libs: ../../Cards/param_card.dat ../../Cards/run_card.dat
	cd ../../Source; make
	touch $@

ifneq (,$(wildcard fbridge.inc))
$(SYCL_BUILDDIR)/.sycllibs:
	$(MAKE) -f $(SYCL_MAKEFILE) commonlib bridgelib_static
	touch $@
endif

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override LIBFLAGSRPATH =
else ifeq ($(USEBUILDDIR),1)
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/../$(LIBDIR)/$(SYCL_BUILDDIR)'
else
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/$(LIBDIR)'
endif


#Build SYCL executable
$(SYCL_BUILDDIR)/$(PROG)_sycl: $(PROCESS) $(DSIG_sycl) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(SYCL_BUILDDIR)/.sycllibs
	$(CXX) $(SYCLFLAGS) -o $(SYCL_BUILDDIR)/$(PROG)_sycl $(PROCESS) $(DSIG_sycl) auto_dsig.o $(MATRIX) $(LIBDIR)/$(SYCL_BUILDDIR)/lib$(SYCL_SYCLLIB).a $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) counters.o -lgomp -lgfortran -L$(LIBDIR)/$(SYCL_BUILDDIR) -l$(SYCL_COMMONLIB) $(LIBFLAGSRPATH) -lstdc++fs
	ln -s $(SYCL_BUILDDIR)/$(PROG)_sycl $(PROG)_sycl

counters.o: counters.cc timer.h
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

$(LIBDIR)libgammaUPC.$(libext):
	cd ../../Source/PDF/gammaUPC; make
endif

# Add source so that the compiler finds the DiscreteSampler module.
$(MATRIX): %.o: %.f
	$(FC) $(FFLAGS) $(MATRIX_FLAG) -c $< -I../../Source/ -I../../Source/PDF/gammaUPC -fopenmp
%.o: %.f
	$(FC) $(FFLAGS) -c $< -I../../Source/ -I../../Source/PDF/gammaUPC -fopenmp -o $@
%_sycl.o: %.f
	$(FC) $(FFLAGS) -c -DMG5AMC_MEEXPORTER_SYCL $< -I../../Source/ -fopenmp -o $@

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
driver_sycl.o: .libs
$(MATRIX): .libs
genps.o: .libs

# Clean

clean:
ifeq (,$(wildcard fbridge.inc))
	$(RM) *.o gensym $(PROG) $(PROG)_forhel
else
	$(RM) *.o gensym $(PROG) $(PROG)_forhel $(SYCL_BUILDDIR)/$(PROG)_sycl $(PROG)_sycl
endif

cleanall:
	make clean
	make -C ../../Source clean
	rm -rf $(LIBDIR)libbias.$(libext)
ifneq (,$(wildcard fbridge.inc))
	$(MAKE) -f $(SYCL_MAKEFILE) cleanall
	rm -f $(SYCL_BUILDDIR)/.sycllibs
endif
	rm -f ../../Source/*.mod ../../Source/*/*.mod
	rm -f .libs
