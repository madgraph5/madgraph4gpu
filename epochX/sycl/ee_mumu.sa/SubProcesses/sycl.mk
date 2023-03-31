#=== Determine the name of this makefile (https://ftp.gnu.org/old-gnu/Manuals/make-3.80/html_node/make_17.html)
#=== NB: different names (e.g. cudacpp.mk and cudacpp_src.mk) are used in the Subprocess and src directories

SYCL_MAKEFILE = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
SYCL_SRC_MAKEFILE = sycl_src.mk

#-------------------------------------------------------------------------------

#=== Use bash in the Makefile (https://www.gnu.org/software/make/manual/html_node/Choosing-the-Shell.html)

SHELL = /bin/bash

#-------------------------------------------------------------------------------

#=== Detect O/S and architecture (assuming uname is available, https://en.wikipedia.org/wiki/Uname)

# Detect O/S kernel (Linux, Darwin...)
UNAME_S := $(shell uname -s)
###$(info UNAME_S='$(UNAME_S)')

# Detect architecture (x86_64, ppc64le...)
UNAME_P := $(shell uname -p)
###$(info UNAME_P='$(UNAME_P)')

#-------------------------------------------------------------------------------

#=== Configure common compiler flags for C++ and SYCL

INCFLAGS = -I.
OPTFLAGS = -O3 -march=native

# Dependency on src directory
MG5AMC_COMMONLIB = mg5amc_common
LIBFLAGS = -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)
INCFLAGS += -I../../src

# Dependency on tools directory
TOOLSDIR = ../../../../../tools
INCFLAGS += -I$(TOOLSDIR)

#-------------------------------------------------------------------------------

#=== Configure the C++ compiler

CXXFLAGS = $(OPTFLAGS) -std=c++20 $(INCFLAGS) -Wall -Wshadow -Wextra
CXXFLAGS+= -ffast-math # see issue #117
ifndef SYCLFLAGS
  $(error SYCLFLAGS not set)
endif

# Note: AR, CXX and FC are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

#-------------------------------------------------------------------------------

#=== Configure ccache for SYCL builds

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for FPTYPE, HELINL, HRDCOD, NTPBMAX

# Set the default FPTYPE (floating point type) choice
ifeq ($(FPTYPE),)
  override FPTYPE = d
endif

# Set the default HELINL (inline helicities?) choice
ifeq ($(HELINL),)
  override HELINL = 0
endif

# Set the default HELINL (inline helicities?) choice
ifeq ($(HRDCOD),)
  override HRDCOD = 1
endif

# Set the default NTPBMAX (maximum threads per block) choice
ifeq ($(NTPBMAX),)
  override NTPBMAX = 1024
endif

# Export FPTYPE, HELINL, HRDCOD, NTPBMAX so that it is not necessary to pass them to the src Makefile too
export FPTYPE
export HELINL
export HRDCOD
export NTPBMAX

#-------------------------------------------------------------------------------

#=== Set the SYCL/C++ compiler flags appropriate to user-defined choices of FPTYPE, HELINL, HRDCOD, NTPBMAX

# Add option to enable CI profiler use
$(info ENABLE_CI_PROFILER=$(ENABLE_CI_PROFILER))
ifeq ($(ENABLE_CI_PROFILER),1)
  CXXFLAGS += --gcc-toolchain="/cvmfs/sft.cern.ch/lcg/releases/gcc/11.3.0-ad0f5/x86_64-centos8"

  # Sets the device ID to the GPU (oneAPI Toolkit) when running check/cmpFcheck in the GitHub CI
  ENABLE_DEVICE_ID = "--device_id=2"
endif

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
$(info FPTYPE=$(FPTYPE))
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'f' and 'd' are supported)
endif

# Set the build flags appropriate to each HELINL choice (example: "make HELINL=1")
$(info HELINL=$(HELINL))
ifeq ($(HELINL),1)
  CXXFLAGS += -DMGONGPU_INLINE_HELAMPS
else ifneq ($(HELINL),0)
  $(error Unknown HELINL='$(HELINL)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each HRDCOD choice (example: "make HRDCOD=1")
$(info HRDCOD=$(HRDCOD))
ifeq ($(HRDCOD),1)
  CXXFLAGS += -DMGONGPU_HARDCODE_PARAM
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each NTPBMAX choice (example: "make NTPBMAX=1024")
$(info NTPBMAX=$(NTPBMAX))
CXXFLAGS += -DMGONGPU_NTPBMAX=$(NTPBMAX)

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
override DIRTAG = $(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
override TAG = $(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
  override LIBDIR = ../../lib/$(BUILDDIR)
  override LIBDIRRPATH = '$$ORIGIN/../$(LIBDIR)'
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is set = 1))
else
  override BUILDDIR = .
  override LIBDIR = ../../lib
  override LIBDIRRPATH = '$$ORIGIN/$(LIBDIR)'
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is not set))
endif

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override CXXLIBFLAGSRPATH =
else
  override CXXLIBFLAGSRPATH = -Wl,-rpath,$(LIBDIRRPATH)
endif

# Setting LD_LIBRARY_PATH or DYLD_LIBRARY_PATH in the RUNTIME is no longer necessary (neither on Linux nor on Mac)
override RUNTIME =

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

sycl_main=$(BUILDDIR)/check.exe
fsycl_main=$(BUILDDIR)/fcheck.exe

# First target (default goal)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(sycl_main) $(fsycl_main)

# Target (and build options): debug
MAKEDEBUG=
debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: MAKEDEBUG := debug
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtagsb)" != "" ]; then echo "Cannot build for tag=$(TAG) as old builds exist for other tags:"; echo "  $(oldtagsb)"; echo "Please run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

# Generic target and build rules: objects from C++ compilation
$(BUILDDIR)/%.o : %.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -c $< -o $@

#-------------------------------------------------------------------------------

# Target (and build rules): common (src) library
commonlib : $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc
	$(MAKE) -C ../../src $(MAKEDEBUG) -f $(SYCL_SRC_MAKEFILE)

#-------------------------------------------------------------------------------

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')

MG5AMC_CXXLIB = mg5amc_$(processid_short)_sycl
cxx_objects_lib=$(BUILDDIR)/CPPProcess.o
cxx_objects_exe=

# Target (and build rules): C++ and SYCL shared libraries
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(BUILDDIR)/fbridge.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: cxx_objects_lib += $(BUILDDIR)/fbridge.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib)
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -shared -o $@ $(cxx_objects_lib) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

# Target (and build rules): C++ and SYCL static libraries
$(LIBDIR)/lib$(MG5AMC_CXXLIB).a: $(BUILDDIR)/CPPProcess.o $(BUILDDIR)/fbridge.o
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	ar rvs $@ $^

#-------------------------------------------------------------------------------

# Target (and build rules): Fortran include files
###$(INCDIR)/%.inc : ../%.inc
###	@if [ ! -d $(INCDIR) ]; then echo "mkdir -p $(INCDIR)"; mkdir -p $(INCDIR); fi
###	\cp $< $@

#-------------------------------------------------------------------------------

# Target (and build rules): C++ and SYCL standalone executables
$(sycl_main): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
###$(sycl_main): $(BUILDDIR)/check_sa.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe)
$(sycl_main): $(LIBDIR)/lib$(MG5AMC_CXXLIB).a $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -o $@ check_sa.cc $(LIBDIR)/lib$(MG5AMC_CXXLIB).a -pthread $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB) -lstdc++fs

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from Fortran compilation
$(BUILDDIR)/%.o : %.f *.inc
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(FC) -I. -c $< -o $@

# Generic target and build rules: objects from Fortran compilation
###$(BUILDDIR)/%.o : %.f *.inc
###	@if [ ! -d $(INCDIR) ]; then echo "mkdir -p $(INCDIR)"; mkdir -p $(INCDIR); fi
###	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
###	$(FC) -I. -I$(INCDIR) -c $< -o $@

# Target (and build rules): Fortran standalone executables
###$(BUILDDIR)/fcheck_sa.o : $(INCDIR)/fbridge.inc

#ifeq ($(UNAME_S),Darwin)
#$(fsycl_main): LIBFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
#endif
#$(fsycl_main): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
#$(fsycl_main): $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).a $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so
#	$(FC) -fsycl -o $@ $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).a $(LIBFLAGS) -lstdc++ -lsycl -L$(LIBDIR) -l$(MG5AMC_COMMONLIB) -lstdc++fs

$(fsycl_main): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(fsycl_main): $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).a $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -o $@ $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).a $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_COMMONLIB) -lstdc++fs

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
ifeq ($(USEBUILDDIR),1)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
	rm -f $(LIBDIR)/lib$(MG5AMC_CXXLIB).so
endif
	$(MAKE) -C ../../src clean -f $(SYCL_SRC_MAKEFILE)
###	rm -rf $(INCDIR)

cleanall:
	@echo
	$(MAKE) USEBUILDDIR=0 clean -f $(SYCL_MAKEFILE)
	@echo
	$(MAKE) USEBUILDDIR=0 -C ../../src cleanall -f $(SYCL_SRC_MAKEFILE)
	rm -rf build.*

# Target: clean the builds as well as the googletest installation
distclean: cleanall
	$(MAKE) -C $(TOOLSDIR) clean
	$(MAKE) -C $(TESTDIR) clean

#-------------------------------------------------------------------------------

# Target: show system and compiler information
info:
	@echo ""
	@uname -spn # e.g. Linux nodename.cern.ch x86_64
ifeq ($(UNAME_S),Darwin)
	@sysctl -a | grep -i brand
	@sysctl -a | grep machdep.cpu | grep features || true
	@sysctl -a | grep hw.physicalcpu:
	@sysctl -a | grep hw.logicalcpu:
else
	@cat /proc/cpuinfo | grep "model name" | sort -u
	@cat /proc/cpuinfo | grep "flags" | sort -u
	@cat /proc/cpuinfo | grep "cpu cores" | sort -u
	@cat /proc/cpuinfo | grep "physical id" | sort -u
endif
	@echo ""
ifneq ($(shell which nvidia-smi 2>/dev/null),)
	nvidia-smi -L
	@echo ""
endif
	@echo USECCACHE=$(USECCACHE)
ifeq ($(USECCACHE),1)
	ccache --version | head -1
endif
	@echo ""
	@echo NVCC=$(NVCC)
ifneq ($(NVCC),)
	$(NVCC) --version
endif
	@echo ""
	@echo CXX=$(CXX)
ifneq ($(shell $(CXX) --version | grep ^clang),)
	@echo $(CXX) -v
	@$(CXX) -v |& egrep -v '(Found|multilib)'
	@readelf -p .comment `$(CXX) -print-libgcc-file-name` |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print "GCC toolchain:",$$5}'
else
	$(CXX) --version
endif
	@echo ""
	@echo FC=$(FC)
	$(FC) --version

#-------------------------------------------------------------------------------

# Target: check (run the C++ test executable)
# [NB THIS IS WHAT IS USED IN THE GITHUB CI!]
check: cmpFcheck

# Target: cmpFcheck (compare ME results from the C++ and Fortran with C++ MEs standalone executables, with a small number of events)
cmpFcheck: all.$(TAG)
	@echo
	@echo "$(BUILDDIR)/check.exe -p 2 32 2 ${ENABLE_DEVICE_ID}"
	@echo "$(BUILDDIR)/fcheck.exe 2 32 2 ${ENABLE_DEVICE_ID}"
	@me1=$(shell $(RUNTIME) $(BUILDDIR)/check.exe -p 2 32 2 ${ENABLE_DEVICE_ID} | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(BUILDDIR)/fcheck.exe 2 32 2 ${ENABLE_DEVICE_ID} | grep Average | awk '{print $$4}'); echo "Avg ME (C++/C++)    = $${me1}"; echo "Avg ME (F77/C++)    = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/C++) returned NaN"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/C++) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi

#-------------------------------------------------------------------------------
