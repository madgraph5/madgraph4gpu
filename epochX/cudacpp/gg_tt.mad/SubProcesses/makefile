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
PLUGIN_MAKEFILE=Makefile
# NB1 Using ":=" below instead of "=" is much faster (it only runs the subprocess once instead of many times)
# NB2 Do not add a comment inlined "PLUGIN_BUILDDIR=$(shell ...) # comment" as otherwise a trailing space is included...
# NB3 The variables relevant to the plugin Makefile must be explicitly passed to $(shell...)
PLUGIN_MAKEENV:=$(shell echo '$(.VARIABLES)' | tr " " "\n" | egrep "(USEBUILDDIR|AVX|FPTYPE|HELINL|HRDCOD)")
###$(info PLUGIN_MAKEENV=$(PLUGIN_MAKEENV))
###$(info $(foreach v,$(PLUGIN_MAKEENV),$(v)="$($(v))"))
PLUGIN_BUILDDIR:=$(shell $(MAKE) $(foreach v,$(PLUGIN_MAKEENV),$(v)="$($(v))") -f $(PLUGIN_MAKEFILE) -pn | awk '/Building/{print $$3}' | sed s/BUILDDIR=//)
###$(info PLUGIN_BUILDDIR='$(PLUGIN_BUILDDIR)')
PLUGIN_COMMONLIB=mg5amc_common
PLUGIN_CXXLIB=mg5amc_$(processid_short)_cpp
PLUGIN_CULIB=mg5amc_$(processid_short)_cuda

LIBS = $(LIBDIR)libbias.$(libext) $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libdsample.$(libext) $(LIBDIR)libgeneric.$(libext) $(LIBDIR)libpdf.$(libext) $(LIBDIR)libmodel.$(libext) $(LIBDIR)libcernlib.$(libext) $(MADLOOP_LIB) $(LOOP_LIBS)

# Source files

MATRIX_HEL = $(patsubst %.f,%.o,$(wildcard matrix*_orig.f))
MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*_optim.f))
ifeq ($(strip $(MATRIX_HEL)),)
        MATRIX = $(patsubst %.f,%.o,$(wildcard matrix*.f))
endif

PROCESS= myamp.o genps.o unwgt.o setcuts.o get_color.o \
         cuts.o cluster.o reweight.o initcluster.o addmothers.o setscales.o \
	 idenparts.o dummy_fct.o

DSIG=driver.o $(patsubst %.f,%.o,$(wildcard auto_dsig*.f))
DSIG_cudacpp=driver_cudacpp.o $(patsubst %.f,%_cudacpp.o,$(wildcard auto_dsig*.f))

SYMMETRY = symmetry.o idenparts.o 

# Binaries

ifeq (,$(wildcard fbridge.inc))
all: $(PROG)
else
all: $(PROG) $(PLUGIN_BUILDDIR)/c$(PROG)_cudacpp $(PLUGIN_BUILDDIR)/g$(PROG)_cudacpp
endif

$(PROG): .libs $(PROCESS) $(DSIG) auto_dsig.o $(LIBS) $(MATRIX) counters.o
	$(FC) -o $(PROG) $(PROCESS) $(DSIG) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o

$(LIBS): .libs

.libs: ../../Cards/param_card.dat ../../Cards/run_card.dat
	cd ../../Source; make
	touch $@

ifneq (,$(wildcard fbridge.inc))
$(PLUGIN_BUILDDIR)/.pluginlibs:
	$(MAKE) -f $(PLUGIN_MAKEFILE)
	touch $@
endif

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override LIBFLAGSRPATH =
else
  override LIBFLAGSRPATH = -Wl,-rpath,'$$ORIGIN/$(LIBDIR)'
endif

$(PLUGIN_BUILDDIR)/c$(PROG)_cudacpp: $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(PLUGIN_BUILDDIR)/.pluginlibs
	$(FC) -o $(PLUGIN_BUILDDIR)/c$(PROG)_cudacpp $(PROCESS) $(DSIG_cudacpp) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o -L$(LIBDIR)/$(PLUGIN_BUILDDIR) -l$(PLUGIN_COMMONLIB) -l$(PLUGIN_CXXLIB) $(LIBFLAGSRPATH)

$(PLUGIN_BUILDDIR)/g$(PROG)_cudacpp: $(PROCESS) $(DSIG_cudacpp) auto_dsig.o $(LIBS) $(MATRIX) counters.o $(PLUGIN_BUILDDIR)/.pluginlibs
	$(FC) -o $(PLUGIN_BUILDDIR)/g$(PROG)_cudacpp $(PROCESS) $(DSIG_cudacpp) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o -L$(LIBDIR)/$(PLUGIN_BUILDDIR) -l$(PLUGIN_COMMONLIB) -l$(PLUGIN_CULIB) $(LIBFLAGSRPATH)

counters.o: counters.cpp timer.h
	$(CXX) -std=c++11 -Wall -Wshadow -Wextra -c $< -o $@

$(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp

gensym: $(SYMMETRY) configs.inc $(LIBDIR)libmodel.$(libext) $(LIBDIR)libgeneric.$(libext)
	$(FC) -o gensym $(SYMMETRY) -L$(LIBDIR) -lmodel -lgeneric -lpdf $(llhapdf) $(LDFLAGS)

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
%_cudacpp.o: %.f
	$(FC) $(FFLAGS) -c -DMG5AMC_MEEXPORTER_CUDACPP $< -I../../Source/ -fopenmp -o $@

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

$(DSIG): .libs
$(DSIG_cudacpp) : .libs
$(MATRIX): .libs
genps.o: .libs

# Plugin avxall targets

ifneq (,$(wildcard fbridge.inc))

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
avxall: avxnone avxsse4
else ifeq ($(UNAME_P),arm)
avxall: avxnone avxsse4
else
avxall: avxnone avxsse4 avxavx2 avx512y avx512z
endif

avxnone: $(PROG) $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=none

avxsse4: $(PROG) $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=sse4

avxavx2: $(PROG) $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=avx2

avx512y: $(PROG) $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512y

avx512z: $(PROG) $(DSIG_cudacpp)
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512z

endif

# Clean

clean:
ifeq (,$(wildcard fbridge.inc))
	$(RM) *.o gensym $(PROG) $(PROG)_forhel
else
	$(RM) *.o gensym $(PROG) $(PROG)_forhel $(PLUGIN_BUILDDIR)/*$(PROG)_cudacpp
endif

cleanall:
	make clean
	make -C ../../Source clean
	rm -rf $(LIBDIR)libbias.$(libext)
ifneq (,$(wildcard fbridge.inc))
	$(MAKE) -f $(PLUGIN_MAKEFILE) cleanall
	rm -f $(PLUGIN_BUILDDIR)/.pluginlibs
endif
	rm -f ../../Source/*.mod ../../Source/*/*.mod
	rm -f .libs
