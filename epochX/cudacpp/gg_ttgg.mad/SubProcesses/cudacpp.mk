# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: S. Hageboeck, J. Teig, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

# This makefile extends the Fortran makefile called "makefile"

CUDACPP_SRC_MAKEFILE = cudacpp_src.mk

# Self-invocation with adapted flags:
cppnative: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=native AVXFLAGS="-march=native" cppbuild
cppnone: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=none AVXFLAGS= cppbuild
cppsse4: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=sse4 AVXFLAGS=-march=nehalem cppbuild
cppavx2: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=avx2 AVXFLAGS=-march=haswell cppbuild
cppavx512y: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=512y AVXFLAGS="-march=skylake-avx512 -mprefer-vector-width=256" cppbuild
cppavx512z: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=512z AVXFLAGS="-march=skylake-avx512 -DMGONGPU_PVW512" cppbuild
cuda: $(SOURCEDIR_GUARD) $(PROCESS)
	$(MAKE) AVX=cuda cudabuild

#-------------------------------------------------------------------------------

#=== Configure common compiler flags for C++ and CUDA
# NB: The base flags are defined in the fortran "makefile"

# Include directories
INCFLAGS = -I. -I../../src

MG_CXXFLAGS  += $(INCFLAGS)
MG_NVCCFLAGS += $(INCFLAGS)

# Dependency on src directory
MG5AMC_COMMONLIB  = mg5amc_common

# Compiler-specific googletest build directory (#125 and #738)
ifneq ($(shell $(CXX) --version | grep '^Intel(R) oneAPI DPC++/C++ Compiler'),)
override CXXNAME = icpx$(shell $(CXX) --version | head -1 | cut -d' ' -f5)
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
override CXXNAME = clang$(shell $(CXX) --version | head -1 | cut -d' ' -f3)
else ifneq ($(shell $(CXX) --version | grep '^g++ (GCC)'),)
override CXXNAME = gcc$(shell $(CXX) --version | head -1 | cut -d' ' -f3)
else
override CXXNAME = unknown
endif
###$(info CXXNAME=$(CXXNAME))
override CXXNAMESUFFIX = _$(CXXNAME)
export CXXNAMESUFFIX

# Dependency on test directory
# Within the madgraph4gpu git repo: by default use a common gtest installation in <topdir>/test (optionally use an external or local gtest)
# Outside the madgraph4gpu git repo: by default do not build the tests (optionally use an external or local gtest)
###GTEST_ROOT = /cvmfs/sft.cern.ch/lcg/releases/gtest/1.11.0-21e8c/x86_64-centos8-gcc11-opt/# example of an external gtest installation
###LOCALGTEST = yes# comment this out (or use make LOCALGTEST=yes) to build tests using a local gtest installation
TESTDIRCOMMON = ../../../../../test
TESTDIRLOCAL = ../../test
ifneq ($(wildcard $(GTEST_ROOT)),)
TESTDIR =
else ifneq ($(LOCALGTEST),)
TESTDIR=$(TESTDIRLOCAL)
GTEST_ROOT = $(TESTDIR)/googletest/install$(CXXNAMESUFFIX)
else ifneq ($(wildcard ../../../../../epochX/cudacpp/CODEGEN),)
TESTDIR = $(TESTDIRCOMMON)
GTEST_ROOT = $(TESTDIR)/googletest/install$(CXXNAMESUFFIX)
else
TESTDIR =
endif
ifneq ($(GTEST_ROOT),)
GTESTLIBDIR = $(GTEST_ROOT)/lib64/
GTESTLIBS = $(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a
GTESTINC = -I$(GTEST_ROOT)/include
else
GTESTLIBDIR =
GTESTLIBS =
GTESTINC =
endif
###$(info GTEST_ROOT = $(GTEST_ROOT))
###$(info LOCALGTEST = $(LOCALGTEST))
###$(info TESTDIR = $(TESTDIR))

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for OMPFLAGS, AVX, FPTYPE, HELINL, HRDCOD, RNDGEN

# Set the default OMPFLAGS choice
OMPFLAGS ?= -fopenmp
ifeq ($(UNAME_S),Darwin) # OM for Mac (any compiler)
override OMPFLAGS = # AV disable OpenMP MT on Apple clang (builds fail in the CI #578)
endif

# Export here, so sub makes don't fall back to the defaults:
export OMPFLAGS

MG_CXXFLAGS += $(OMPFLAGS)

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the build directory)
DIRTAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)
CUDACPP_BUILDDIR = build.$(DIRTAG)
CUDACPP_LIBDIR := ../../lib/$(CUDACPP_BUILDDIR)
LIBDIRRPATH := '$$ORIGIN:$$ORIGIN/../$(CUDACPP_LIBDIR)'
ifneq ($(AVX),)
  $(info Building CUDACPP in CUDACPP_BUILDDIR=$(CUDACPP_BUILDDIR). Libs in $(CUDACPP_LIBDIR))
endif

# On Linux, set rpath to CUDACPP_LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables or shared libraries ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in CUDACPP_LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override CXXLIBFLAGSRPATH =
  override CULIBFLAGSRPATH =
else
  # RPATH to cuda/cpp libs when linking executables
  override CXXLIBFLAGSRPATH = -Wl,-rpath,$(LIBDIRRPATH)
  override CULIBFLAGSRPATH = -Xlinker -rpath,$(LIBDIRRPATH)
endif

# Setting LD_LIBRARY_PATH or DYLD_LIBRARY_PATH in the RUNTIME is no longer necessary (neither on Linux nor on Mac)
override RUNTIME =

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

cxx_main=$(CUDACPP_BUILDDIR)/check.exe
fcxx_main=$(CUDACPP_BUILDDIR)/fcheck.exe

cu_main=$(CUDACPP_BUILDDIR)/gcheck.exe
fcu_main=$(CUDACPP_BUILDDIR)/fgcheck.exe

ifneq ($(GTESTLIBS),)
testmain=$(CUDACPP_BUILDDIR)/runTest.exe
cutestmain=$(CUDACPP_BUILDDIR)/runTest_cuda.exe
endif

cppbuild: $(CUDACPP_BUILDDIR)/$(PROG)_cpp $(cxx_main) $(fcxx_main) $(testmain)
cudabuild: $(CUDACPP_BUILDDIR)/$(PROG)_cuda $(cu_main) $(fcu_main) $(cutestmain)

# Generic target and build rules: objects from CUDA compilation
$(CUDACPP_BUILDDIR)/%.o : %.cu *.h ../../src/*.h
	@mkdir -p $(CUDACPP_BUILDDIR)
	$(NVCC) $(MG_NVCCFLAGS) $(NVCCFLAGS) -c $< -o $@

$(CUDACPP_BUILDDIR)/%_cu.o : %.cc *.h ../../src/*.h
	@mkdir -p $(CUDACPP_BUILDDIR)
	$(NVCC) $(MG_NVCCFLAGS) $(NVCCFLAGS) -c -x cu $< -o $@

# Generic target and build rules: objects from C++ compilation
# (NB do not include CUINC here! add it only for NVTX or curand #679)
$(CUDACPP_BUILDDIR)/%.o : %.cc *.h ../../src/*.h
	@mkdir -p $(CUDACPP_BUILDDIR)
	$(CXX) $(MG_CXXFLAGS) $(CXXFLAGS) -c $< -o $@

# Apply special build flags only to CrossSectionKernel.cc and gCrossSectionKernel.cu (no fast math, see #117 and #516)
ifeq ($(shell $(CXX) --version | grep ^nvc++),)
$(CUDACPP_BUILDDIR)/CrossSectionKernels.o: CXXFLAGS += -fno-fast-math
ifneq ($(NVCC),)
$(CUDACPP_BUILDDIR)/gCrossSectionKernels.o: NVCCFLAGS += -Xcompiler -fno-fast-math
endif
endif

# Apply special build flags only to check_sa.o and gcheck_sa.o (NVTX in timermap.h, #679)
$(CUDACPP_BUILDDIR)/check_sa.o: MG_CXXFLAGS += $(USE_NVTX) $(CUINC)
$(CUDACPP_BUILDDIR)/gcheck_sa.o: MG_CXXFLAGS += $(USE_NVTX) $(CUINC)

# Apply special build flags only to check_sa and CurandRandomNumberKernel (curand headers, #679)
$(CUDACPP_BUILDDIR)/check_sa.o: MG_CXXFLAGS += $(CXXFLAGSCURAND)
$(CUDACPP_BUILDDIR)/gcheck_sa.o: MG_NVCCFLAGS += $(CXXFLAGSCURAND)
$(CUDACPP_BUILDDIR)/CurandRandomNumberKernel.o: MG_CXXFLAGS += $(CXXFLAGSCURAND)
$(CUDACPP_BUILDDIR)/gCurandRandomNumberKernel.o: MG_NVCCFLAGS += $(CXXFLAGSCURAND)


# Avoid "warning: builtin __has_trivial_... is deprecated; use __is_trivially_... instead" in nvcc with icx2023 (#592)
ifneq ($(shell $(CXX) --version | egrep '^(Intel)'),)
ifneq ($(NVCC),)
MG_NVCCFLAGS += -Xcompiler -Wno-deprecated-builtins
endif
endif

#### Apply special build flags only to CPPProcess.cc (-flto)
###$(BUILDDIR)/CPPProcess.o: CXXFLAGS += -flto

#-------------------------------------------------------------------------------

$(CUDACPP_LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc
	$(MAKE) AVX=$(AVX) AVXFLAGS="$(AVXFLAGS)" -C ../../src -f $(CUDACPP_SRC_MAKEFILE)

#-------------------------------------------------------------------------------

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')
###$(info processid_short=$(processid_short))

MG5AMC_CXXLIB = mg5amc_$(processid_short)_cpp
cxx_objects_lib=$(CUDACPP_BUILDDIR)/CPPProcess.o $(CUDACPP_BUILDDIR)/MatrixElementKernels.o $(CUDACPP_BUILDDIR)/BridgeKernels.o $(CUDACPP_BUILDDIR)/CrossSectionKernels.o
cxx_objects_exe=$(CUDACPP_BUILDDIR)/CommonRandomNumberKernel.o $(CUDACPP_BUILDDIR)/RamboSamplingKernels.o

MG5AMC_CULIB = mg5amc_$(processid_short)_cuda
cu_objects_lib=$(CUDACPP_BUILDDIR)/gCPPProcess.o $(CUDACPP_BUILDDIR)/gMatrixElementKernels.o $(CUDACPP_BUILDDIR)/gBridgeKernels.o $(CUDACPP_BUILDDIR)/gCrossSectionKernels.o
cu_objects_exe=$(CUDACPP_BUILDDIR)/gCommonRandomNumberKernel.o $(CUDACPP_BUILDDIR)/gRamboSamplingKernels.o

# Target (and build rules): C++ and CUDA shared libraries
$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so: $(CUDACPP_BUILDDIR)/fbridge.o
$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so: cxx_objects_lib += $(CUDACPP_BUILDDIR)/fbridge.o
$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so: $(CUDACPP_LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib)
	$(CXX) -shared -o $@ $(cxx_objects_lib) $(CXXLIBFLAGSRPATH) -L$(CUDACPP_LIBDIR) -l$(MG5AMC_COMMONLIB) $(MG_LDFLAGS) $(LDFLAGS)

$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so: $(CUDACPP_BUILDDIR)/fbridge_cu.o
$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so: cu_objects_lib += $(CUDACPP_BUILDDIR)/fbridge_cu.o
$(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so: $(CUDACPP_LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_objects_lib)
	$(NVCC) --shared -o $@ $(cu_objects_lib) $(CULIBFLAGSRPATH) -L$(CUDACPP_LIBDIR) -l$(MG5AMC_COMMONLIB)

#-------------------------------------------------------------------------------

# Target (and build rules): C++ and CUDA standalone executables

$(cxx_main): MG_LDFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cxx_main): MG_LDFLAGS += -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CXXLIB) # Process-specific library
$(cxx_main): $(CUDACPP_BUILDDIR)/check_sa.o $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe) $(CUDACPP_BUILDDIR)/CurandRandomNumberKernel.o
	$(CXX) -o $@ $(CUDACPP_BUILDDIR)/check_sa.o $(OMPFLAGS) -ldl -pthread $(cxx_objects_exe) $(CUDACPP_BUILDDIR)/CurandRandomNumberKernel.o $(MG_LDFLAGS) $(LDFLAGS)

ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(cu_main): MG_LDFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(cu_main): MG_LDFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(cu_main): MG_LDFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif
$(cu_main): MG_LDFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cu_main): MG_LDFLAGS += -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CULIB) # Process-specific library
$(cu_main): $(CUDACPP_BUILDDIR)/gcheck_sa.o $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so $(cu_objects_exe) $(CUDACPP_BUILDDIR)/gCurandRandomNumberKernel.o
	$(NVCC) -o $@ $(CUDACPP_BUILDDIR)/gcheck_sa.o $(CUARCHFLAGS) $(cu_objects_exe) $(CUDACPP_BUILDDIR)/gCurandRandomNumberKernel.o $(MG_LDFLAGS) $(LDFLAGS)

#-------------------------------------------------------------------------------
# Check executables:

ifeq ($(UNAME_S),Darwin)
$(fcxx_main): MG_LDFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(fcxx_main): MG_LDFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(fcxx_main): MG_LDFLAGS += -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CXXLIB) # Process-specific library
$(fcxx_main): $(CUDACPP_BUILDDIR)/fcheck_sa.o $(CUDACPP_BUILDDIR)/fsampler.o $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe)
	$(CXX) -o $@ $(CUDACPP_BUILDDIR)/fcheck_sa.o $(cxx_objects_exe) $(OMPFLAGS) $(CUDACPP_BUILDDIR)/fsampler.o -lgfortran -L$(CUDACPP_LIBDIR) $(MG_LDFLAGS) $(LDFLAGS)

ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(fcu_main): MG_LDFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(fcu_main): MG_LDFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
endif
ifeq ($(UNAME_S),Darwin)
$(fcu_main): MG_LDFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(fcu_main): MG_LDFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(fcu_main): MG_LDFLAGS += -L$(CUDACPP_BUILDDIR) -l$(MG5AMC_CULIB) # Process-specific library
$(fcu_main): $(CUDACPP_BUILDDIR)/fcheck_sa.o $(CUDACPP_BUILDDIR)/fsampler_cu.o $(CUDACPP_BUILDDIR)/lib$(MG5AMC_CULIB).so $(cu_objects_exe)
	$(NVCC) -o $@ $(CUDACPP_BUILDDIR)/fcheck_sa.o $(CUDACPP_BUILDDIR)/fsampler_cu.o $(cu_objects_exe)  -lgfortran $(MG_LDFLAGS) $(LDFLAGS)

#-------------------------------------------------------------------------------

# Target (and build rules): test objects and test executable

$(testmain) $(cutestmain): $(GTESTLIBS)
$(testmain) $(cutestmain): INCFLAGS += $(GTESTINC)
$(testmain) $(cutestmain): MG_LDFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main

$(CUDACPP_BUILDDIR)/testxxx.o $(CUDACPP_BUILDDIR)/testxxx_cu.o: $(GTESTLIBS) testxxx_cc_ref.txt
$(testmain): $(CUDACPP_BUILDDIR)/testxxx.o
$(testmain): cxx_objects_exe += $(CUDACPP_BUILDDIR)/testxxx.o # Comment out this line to skip the C++ test of xxx functions
$(cutestmain): $(CUDACPP_BUILDDIR)/testxxx_cu.o
$(cutestmain): cu_objects_exe += $(CUDACPP_BUILDDIR)/testxxx_cu.o # Comment out this line to skip the CUDA test of xxx functions


$(CUDACPP_BUILDDIR)/testmisc.o $(CUDACPP_BUILDDIR)/testmisc_cu.o: $(GTESTLIBS)
$(testmain): $(CUDACPP_BUILDDIR)/testmisc.o
$(testmain): cxx_objects_exe += $(CUDACPP_BUILDDIR)/testmisc.o # Comment out this line to skip the C++ miscellaneous tests
$(cutestmain): $(CUDACPP_BUILDDIR)/testmisc_cu.o
$(cutestmain): cu_objects_exe += $(CUDACPP_BUILDDIR)/testmisc_cu.o # Comment out this line to skip the CUDA miscellaneous tests


$(CUDACPP_BUILDDIR)/runTest.o $(CUDACPP_BUILDDIR)/runTest_cu.o: $(GTESTLIBS)
$(testmain): $(CUDACPP_BUILDDIR)/runTest.o
$(testmain): cxx_objects_exe += $(CUDACPP_BUILDDIR)/runTest.o
$(cutestmain): $(CUDACPP_BUILDDIR)/runTest_cu.o
$(cutestmain): cu_objects_exe  += $(CUDACPP_BUILDDIR)/runTest_cu.o


ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(cutestmain): MG_LDFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(cutestmain): MG_LDFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(cutestmain): MG_LDFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif


ifneq ($(OMPFLAGS),)
ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
$(testmain): MG_LDFLAGS += -liomp5 # see #578 (not '-qopenmp -static-intel' as in https://stackoverflow.com/questions/45909648)
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
$(testmain): MG_LDFLAGS += -L $(shell dirname $(shell $(CXX) -print-file-name=libc++.so)) -lomp # see #604
###else ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
###$(testmain): LIBFLAGS += ???? # OMP is not supported yet by cudacpp for Apple clang (see #578 and #604)
else
$(testmain): MG_LDFLAGS += -lgomp
endif
endif

$(testmain): MG_LDFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(testmain): $(CUDACPP_LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib) $(cxx_objects_exe) $(GTESTLIBS)
	$(CXX) -o $@ $(cxx_objects_lib) $(cxx_objects_exe) -L$(CUDACPP_LIBDIR) -l$(MG5AMC_COMMONLIB) -ldl -pthread $(MG_LDFLAGS) $(LDFLAGS)

$(cutestmain): MG_LDFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cutestmain): $(CUDACPP_LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_objects_lib) $(cu_objects_exe) $(GTESTLIBS)
	$(NVCC) -o $@ $(cu_objects_lib) $(cu_objects_exe) -L$(CUDACPP_LIBDIR) -l$(MG5AMC_COMMONLIB) -ldl -lcuda $(MG_LDFLAGS) $(LDFLAGS)

# Use target gtestlibs to build only googletest
ifneq ($(GTESTLIBS),)
gtestlibs: $(GTESTLIBS)
endif

# Use flock (Linux only, no Mac) to allow 'make -j' if googletest has not yet been downloaded https://stackoverflow.com/a/32666215
$(GTESTLIBS):
ifneq ($(shell which flock 2>/dev/null),)
	flock $(TESTDIR)/.make_test.lock $(MAKE) -C $(TESTDIR)
else
	if [ -d $(TESTDIR) ]; then $(MAKE) -C $(TESTDIR); fi
endif

#-------------------------------------------------------------------------------

# Target: clean the builds as well as the gtest installation(s)
distclean: clean cleansrc
ifneq ($(wildcard $(TESTDIRCOMMON)),)
	$(MAKE) -C $(TESTDIRCOMMON) clean
endif
	$(MAKE) -C $(TESTDIRLOCAL) clean

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

# Target: check/gcheck (run the C++ test executable)
# [NB THIS IS WHAT IS USED IN THE GITHUB CI!]
check: runTest cmpFcheck
gcheck:
	$(MAKE) AVX=cuda runTest cmpFGcheck

# Target: runTest (run the C++ test executable runTest.exe)
ifneq ($(AVX),cuda)
runTest: cppbuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/runTest.exe
else
runTest: cudabuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/runTest_cuda.exe
endif


# Target: runCheck (run the C++ standalone executable check.exe, with a small number of events)
runCheck: cppbuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/check.exe -p 2 32 2

# Target: runGcheck (run the CUDA standalone executable gcheck.exe, with a small number of events)
runGcheck: AVX=cuda
runGcheck:
	$(MAKE) AVX=cuda cudabuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/gcheck.exe -p 2 32 2

# Target: runFcheck (run the Fortran standalone executable - with C++ MEs - fcheck.exe, with a small number of events)
runFcheck: cppbuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/fcheck.exe 2 32 2

# Target: runFGcheck (run the Fortran standalone executable - with CUDA MEs - fgcheck.exe, with a small number of events)
runFGcheck: AVX=cuda
runFGcheck:
	$(MAKE) AVX=cuda cudabuild
	$(RUNTIME) $(CUDACPP_BUILDDIR)/fgcheck.exe 2 32 2

# Target: cmpFcheck (compare ME results from the C++ and Fortran with C++ MEs standalone executables, with a small number of events)
cmpFcheck: cppbuild
	@echo
	@echo "$(CUDACPP_BUILDDIR)/check.exe --common -p 2 32 2"
	@echo "$(CUDACPP_BUILDDIR)/fcheck.exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(CUDACPP_BUILDDIR)/check.exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(CUDACPP_BUILDDIR)/fcheck.exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/C++)    = $${me1}"; echo "Avg ME (F77/C++)    = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/C++) returned NaN"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/C++) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi

# Target: cmpFGcheck (compare ME results from the CUDA and Fortran with CUDA MEs standalone executables, with a small number of events)
cmpFGcheck: AVX=cuda
cmpFGcheck:
	$(MAKE) AVX=cuda cudabuild
	@echo
	@echo "$(CUDACPP_BUILDDIR)/gcheck.exe --common -p 2 32 2"
	@echo "$(CUDACPP_BUILDDIR)/fgcheck.exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(CUDACPP_BUILDDIR)/gcheck.exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(CUDACPP_BUILDDIR)/fgcheck.exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/CUDA)   = $${me1}"; echo "Avg ME (F77/CUDA)   = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/CUDA) crashed"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/CUDA) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi

