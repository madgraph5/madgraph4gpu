# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

#=== Determine the name of this makefile (https://ftp.gnu.org/old-gnu/Manuals/make-3.80/html_node/make_17.html)
#=== NB: use ':=' to ensure that the value of CUDACPP_MAKEFILE is not modified further down after including make_opts
#=== NB: use 'override' to ensure that the value can not be modified from the outside
override CUDACPP_MAKEFILE := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
###$(info CUDACPP_MAKEFILE='$(CUDACPP_MAKEFILE)')

#=== NB: different names (e.g. cudacpp.mk and cudacpp_src.mk) are used in the Subprocess and src directories
override CUDACPP_SRC_MAKEFILE = cudacpp_src.mk

#-------------------------------------------------------------------------------

#=== Use bash in the Makefile (https://www.gnu.org/software/make/manual/html_node/Choosing-the-Shell.html)

SHELL := /bin/bash

#-------------------------------------------------------------------------------

#=== Detect O/S and architecture (assuming uname is available, https://en.wikipedia.org/wiki/Uname)

# Detect O/S kernel (Linux, Darwin...)
UNAME_S := $(shell uname -s)
###$(info UNAME_S='$(UNAME_S)')

# Detect architecture (x86_64, ppc64le...)
UNAME_P := $(shell uname -p)
###$(info UNAME_P='$(UNAME_P)')

#-------------------------------------------------------------------------------

#=== Include the common MG5aMC Makefile options

# OM: this is crucial for MG5aMC flag consistency/documentation
# AV: temporarely comment this out because it breaks cudacpp builds
ifneq ($(wildcard ../../Source/make_opts),)
include ../../Source/make_opts
endif

#-------------------------------------------------------------------------------

#=== Configure common compiler flags for C++ and CUDA

INCFLAGS = -I.
OPTFLAGS = -O3 # this ends up in CUFLAGS too (should it?), cannot add -Ofast or -ffast-math here

# Dependency on src directory
MG5AMC_COMMONLIB = mg5amc_common
LIBFLAGS = -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)
INCFLAGS += -I../../src

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
TESTDIR =
GTESTLIBDIR =
GTESTLIBS =
GTESTINC =

###$(info GTEST_ROOT = $(GTEST_ROOT))
###$(info LOCALGTEST = $(LOCALGTEST))
###$(info TESTDIR = $(TESTDIR))

TESTDIR=
GTESTLIBDIR =
GTESTLIBS =
GTESTINC =
#-------------------------------------------------------------------------------

#=== Configure the C++ compiler

CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) -Wall -Wshadow -Wextra
ifeq ($(shell $(CXX) --version | grep ^nvc++),)
CXXFLAGS += -ffast-math # see issue #117
endif
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
###CXXFLAGS+= -g # FOR DEBUGGING ONLY

# Optionally add debug flags to display the full list of flags (eg on Darwin)
###CXXFLAGS+= -v

# Note: AR, CXX and FC are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

# Add -mmacosx-version-min=11.3 to avoid "ld: warning: object file was built for newer macOS version than being linked"
ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
CXXFLAGS += -mmacosx-version-min=11.3
endif

#-------------------------------------------------------------------------------

#=== Configure the CUDA compiler

# If CXX is not a single word (example "clang++ --gcc-toolchain...") then disable CUDA builds (issue #505)
# This is because it is impossible to pass this to "CUFLAGS += -ccbin <host-compiler>" below
ifneq ($(words $(subst ccache ,,$(CXX))),1) # allow at most "CXX=ccache <host-compiler>" from outside
  $(warning CUDA builds are not supported for multi-word CXX "$(CXX)")
  override CUDA_HOME=disabled
endif

# If CUDA_HOME is not set, try to set it from the location of nvcc
ifndef CUDA_HOME
  CUDA_HOME = $(patsubst %bin/nvcc,%,$(shell which nvcc 2>/dev/null))
  $(warning CUDA_HOME was not set: using "$(CUDA_HOME)")
endif

# Set NVCC as $(CUDA_HOME)/bin/nvcc if it exists
ifneq ($(wildcard $(CUDA_HOME)/bin/nvcc),)
  NVCC = $(CUDA_HOME)/bin/nvcc
  USE_NVTX ?=-DUSE_NVTX
  # See https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
  # See https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
  # Default: use compute capability 70 for V100 (CERN lxbatch, CERN itscrd, Juwels Cluster).
  # Embed device code for 70, and PTX for 70+.
  # Export MADGRAPH_CUDA_ARCHITECTURE (comma-separated list) to use another value or list of values (see #533).
  # Examples: use 60 for P100 (Piz Daint), 80 for A100 (Juwels Booster, NVidia raplab/Curiosity).
  MADGRAPH_CUDA_ARCHITECTURE ?= 70
  ###CUARCHFLAGS = -gencode arch=compute_$(MADGRAPH_CUDA_ARCHITECTURE),code=compute_$(MADGRAPH_CUDA_ARCHITECTURE) -gencode arch=compute_$(MADGRAPH_CUDA_ARCHITECTURE),code=sm_$(MADGRAPH_CUDA_ARCHITECTURE) # Older implementation (AV): go back to this one for multi-GPU support #533
  ###CUARCHFLAGS = --gpu-architecture=compute_$(MADGRAPH_CUDA_ARCHITECTURE) --gpu-code=sm_$(MADGRAPH_CUDA_ARCHITECTURE),compute_$(MADGRAPH_CUDA_ARCHITECTURE) # Newer implementation (SH): cannot use this as-is for multi-GPU support #533
  comma:=,
  CUARCHFLAGS = $(foreach arch,$(subst $(comma), ,$(MADGRAPH_CUDA_ARCHITECTURE)),-gencode arch=compute_$(arch),code=compute_$(arch) -gencode arch=compute_$(arch),code=sm_$(arch))
  CUINC = -I$(CUDA_HOME)/include/
  ifeq ($(RNDGEN),hasNoCurand)
    CURANDLIBFLAGS=
  else
    CURANDLIBFLAGS = -L$(CUDA_HOME)/lib64/ -lcurand # NB: -lcuda is not needed here!
  endif
  CUOPTFLAGS = -lineinfo
  CUFLAGS = $(foreach opt, $(OPTFLAGS), -Xcompiler $(opt)) $(CUOPTFLAGS) $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math
  ###CUFLAGS += -Xcompiler -Wall -Xcompiler -Wextra -Xcompiler -Wshadow
  ###NVCC_VERSION = $(shell $(NVCC) --version | grep 'Cuda compilation tools' | cut -d' ' -f5 | cut -d, -f1)
  CUFLAGS += -std=c++17 # need CUDA >= 11.2 (see #333): this is enforced in mgOnGpuConfig.h
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)
else ifneq ($(origin REQUIRE_CUDA),undefined)
  # If REQUIRE_CUDA is set but no cuda is found, stop here (e.g. for CI tests on GPU #443)
  $(error No cuda installation found (set CUDA_HOME or make nvcc visible in PATH))
else
  # No cuda. Switch cuda compilation off and go to common random numbers in C++
  $(warning CUDA_HOME is not set or is invalid: export CUDA_HOME to compile with cuda)
  override NVCC=
  override USE_NVTX=
  override CUINC=
  override CURANDLIBFLAGS=
endif
export NVCC
export CUFLAGS

# Set the host C++ compiler for nvcc via "-ccbin <host-compiler>"
# (NB issue #505: this must be a single word, "clang++ --gcc-toolchain..." is not supported)
CUFLAGS += -ccbin $(shell which $(subst ccache ,,$(CXX)))

# Allow newer (unsupported) C++ compilers with older versions of CUDA if ALLOW_UNSUPPORTED_COMPILER_IN_CUDA is set (#504)
ifneq ($(origin ALLOW_UNSUPPORTED_COMPILER_IN_CUDA),undefined)
CUFLAGS += -allow-unsupported-compiler
endif

#-------------------------------------------------------------------------------

#=== Configure ccache for C++ and CUDA builds

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif
#ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
#  override AR:=ccache $(AR)
#endif
ifneq ($(NVCC),)
  ifeq ($(USECCACHE)$(shell echo $(NVCC) | grep ccache),1)
    override NVCC:=ccache $(NVCC)
  endif
endif

#-------------------------------------------------------------------------------

#=== Configure PowerPC-specific compiler flags for C++ and CUDA

# PowerPC-specific CXX compiler flags (being reviewed)
ifeq ($(UNAME_P),ppc64le)
  CXXFLAGS+= -mcpu=power9 -mtune=power9 # gains ~2-3% both for none and sse4
  # Throughput references without the extra flags below: none=1.41-1.42E6, sse4=2.15-2.19E6
  ###CXXFLAGS+= -DNO_WARN_X86_INTRINSICS # no change
  ###CXXFLAGS+= -fpeel-loops # no change
  ###CXXFLAGS+= -funroll-loops # gains ~1% for none, loses ~1% for sse4
  ###CXXFLAGS+= -ftree-vectorize # no change
  ###CXXFLAGS+= -flto # would increase to none=4.08-4.12E6, sse4=4.99-5.03E6!
else
  ###CXXFLAGS+= -flto # also on Intel this would increase throughputs by a factor 2 to 4...
  ######CXXFLAGS+= -fno-semantic-interposition # no benefit (neither alone, nor combined with -flto)
endif

# PowerPC-specific CUDA compiler flags (to be reviewed!)
ifeq ($(UNAME_P),ppc64le)
  CUFLAGS+= -Xcompiler -mno-float128
endif

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for OMPFLAGS, AVX, FPTYPE, HELINL, HRDCOD, RNDGEN

# Set the default OMPFLAGS choice
ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
override OMPFLAGS = -fopenmp
###override OMPFLAGS = # disable OpenMP MT on Intel (was ok without nvcc but not ok with nvcc before #578)
else ifneq ($(shell $(CXX) --version | egrep '^(clang)'),)
override OMPFLAGS = -fopenmp
###override OMPFLAGS = # disable OpenMP MT on clang (was not ok without or with nvcc before #578)
###else ifneq ($(shell $(CXX) --version | egrep '^(Apple clang)'),) # AV for Mac (Apple clang compiler)
else ifeq ($(UNAME_S),Darwin) # OM for Mac (any compiler)
override OMPFLAGS = # AV disable OpenMP MT on Apple clang (builds fail in the CI #578)
###override OMPFLAGS = -fopenmp # OM reenable OpenMP MT on Apple clang? (AV Oct 2023: this still fails in the CI)
else
override OMPFLAGS = -fopenmp # enable OpenMP MT by default on all other platforms
###override OMPFLAGS = # disable OpenMP MT on all other platforms (default before #575)
endif

# Set the default AVX (vectorization) choice
ifeq ($(AVX),)
  ifeq ($(UNAME_P),ppc64le)
    ###override AVX = none
    override AVX = sse4
  else ifeq ($(UNAME_P),arm)
    ###override AVX = none
    override AVX = sse4
  else ifeq ($(wildcard /proc/cpuinfo),)
    override AVX = none
    $(warning Using AVX='$(AVX)' because host SIMD features cannot be read from /proc/cpuinfo)
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

# Set the default FPTYPE (floating point type) choice
ifeq ($(FPTYPE),)
  override FPTYPE = d
endif

# Set the default HELINL (inline helicities?) choice
ifeq ($(HELINL),)
  override HELINL = 0
endif

# Set the default HRDCOD (hardcode cIPD physics parameters?) choice
ifeq ($(HRDCOD),)
  override HRDCOD = 0
endif

# Set the default RNDGEN (random number generator) choice
ifeq ($(RNDGEN),)
  ifeq ($(NVCC),)
    override RNDGEN = hasNoCurand
  else ifeq ($(RNDGEN),)
    override RNDGEN = hasCurand
  endif
endif

# Export AVX, FPTYPE, HELINL, HRDCOD, RNDGEN, OMPFLAGS so that it is not necessary to pass them to the src Makefile too
export AVX
export FPTYPE
export HELINL
export HRDCOD
export RNDGEN
export OMPFLAGS

#-------------------------------------------------------------------------------

#=== Set the CUDA/C++ compiler flags appropriate to user-defined choices of AVX, FPTYPE, HELINL, HRDCOD, RNDGEN

# Set the build flags appropriate to OMPFLAGS
$(info OMPFLAGS=$(OMPFLAGS))
CXXFLAGS += $(OMPFLAGS)

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
else ifeq ($(UNAME_P),arm)
  ifeq ($(AVX),sse4)
    override AVXFLAGS = -D__SSE4_2__ # ARM NEON with 128 width (Q/quadword registers)
  else ifneq ($(AVX),none)
    $(error Unknown AVX='$(AVX)': only 'none' and 'sse4' are supported on ARM for the moment)
  endif
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
  ifeq ($(AVX),none)
    override AVXFLAGS = -mno-sse3 # no SIMD
  else ifeq ($(AVX),sse4)
    override AVXFLAGS = -mno-avx # SSE4.2 with 128 width (xmm registers)
  else ifeq ($(AVX),avx2)
    override AVXFLAGS = -march=haswell # AVX2 with 256 width (ymm registers) [DEFAULT for clang]
  else ifeq ($(AVX),512y)
    override AVXFLAGS = -march=skylake -mprefer-vector-width=256 # AVX512 with 256 width (ymm registers) [DEFAULT for gcc]
  else ifeq ($(AVX),512z)
    override AVXFLAGS = -march=skylake -DMGONGPU_PVW512 # AVX512 with 512 width (zmm registers)
  else
    $(error Unknown AVX='$(AVX)': only 'none', 'sse4', 'avx2', '512y' and '512z' are supported)
  endif
else
  ifeq ($(AVX),none)
    override AVXFLAGS = -march=x86-64 # no SIMD (see #588)
  else ifeq ($(AVX),sse4)
    override AVXFLAGS = -march=nehalem # SSE4.2 with 128 width (xmm registers)
  else ifeq ($(AVX),avx2)
    override AVXFLAGS = -march=haswell # AVX2 with 256 width (ymm registers) [DEFAULT for clang]
  else ifeq ($(AVX),512y)
    override AVXFLAGS = -march=skylake-avx512 -mprefer-vector-width=256 # AVX512 with 256 width (ymm registers) [DEFAULT for gcc]
  else ifeq ($(AVX),512z)
    override AVXFLAGS = -march=skylake-avx512 -DMGONGPU_PVW512 # AVX512 with 512 width (zmm registers)
  else
    $(error Unknown AVX='$(AVX)': only 'none', 'sse4', 'avx2', '512y' and '512z' are supported)
  endif
endif
# For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations?
CXXFLAGS+= $(AVXFLAGS)

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
$(info FPTYPE=$(FPTYPE))
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_DOUBLE
  CUFLAGS  += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT -DMGONGPU_FPTYPE2_FLOAT
  CUFLAGS  += -DMGONGPU_FPTYPE_FLOAT -DMGONGPU_FPTYPE2_FLOAT
else ifeq ($(FPTYPE),m)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_FLOAT
  CUFLAGS  += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'd', 'f' and 'm' are supported)
endif

# Set the build flags appropriate to each HELINL choice (example: "make HELINL=1")
$(info HELINL=$(HELINL))
ifeq ($(HELINL),1)
  CXXFLAGS += -DMGONGPU_INLINE_HELAMPS
  CUFLAGS  += -DMGONGPU_INLINE_HELAMPS
else ifneq ($(HELINL),0)
  $(error Unknown HELINL='$(HELINL)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each HRDCOD choice (example: "make HRDCOD=1")
$(info HRDCOD=$(HRDCOD))
ifeq ($(HRDCOD),1)
  CXXFLAGS += -DMGONGPU_HARDCODE_PARAM
  CUFLAGS  += -DMGONGPU_HARDCODE_PARAM
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each RNDGEN choice (example: "make RNDGEN=hasNoCurand")
$(info RNDGEN=$(RNDGEN))
ifeq ($(RNDGEN),hasNoCurand)
  override CXXFLAGSCURAND = -DMGONGPU_HAS_NO_CURAND
else ifeq ($(RNDGEN),hasCurand)
  override CXXFLAGSCURAND =
else
  $(error Unknown RNDGEN='$(RNDGEN)': only 'hasCurand' and 'hasNoCurand' are supported)
endif

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
override DIRTAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of CUDA and no-CUDA environment builds with different random number generators)
override TAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)_$(RNDGEN)

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
###override INCDIR = ../../include
###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables or shared libraries ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override CXXLIBFLAGSRPATH =
  override CULIBFLAGSRPATH =
  override CXXLIBFLAGSRPATH2 =
  override CULIBFLAGSRPATH2 =
else
  # RPATH to cuda/cpp libs when linking executables
  override CXXLIBFLAGSRPATH = -Wl,-rpath,$(LIBDIRRPATH)
  override CULIBFLAGSRPATH = -Xlinker -rpath,$(LIBDIRRPATH)
  # RPATH to common lib when linking cuda/cpp libs
  override CXXLIBFLAGSRPATH2 = -Wl,-rpath,'$$ORIGIN'
  override CULIBFLAGSRPATH2 = -Xlinker -rpath,'$$ORIGIN'
endif

# Setting LD_LIBRARY_PATH or DYLD_LIBRARY_PATH in the RUNTIME is no longer necessary (neither on Linux nor on Mac)
override RUNTIME =

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

cxx_main=$(BUILDDIR)/check.exe
fcxx_main=$(BUILDDIR)/fcheck.exe

ifneq ($(NVCC),)
cu_main=$(BUILDDIR)/gcheck.exe
fcu_main=$(BUILDDIR)/fgcheck.exe
else
cu_main=
fcu_main=
endif

testmain=$(BUILDDIR)/runTest.exe

ifneq ($(GTESTLIBS),)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_main) $(cxx_main) $(fcu_main) $(fcxx_main) $(testmain)
else
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_main) $(cxx_main) $(fcu_main) $(fcxx_main)
endif

# Target (and build options): debug
MAKEDEBUG=
debug: OPTFLAGS   = -g -O0
debug: CUOPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtagsb)" != "" ]; then echo "Cannot build for tag=$(TAG) as old builds exist for other tags:"; echo "  $(oldtagsb)"; echo "Please run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

# Generic target and build rules: objects from CUDA compilation
ifneq ($(NVCC),)
$(BUILDDIR)/%.o : %.cu *.h ../../src/*.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -Xcompiler -fPIC -c $< -o $@

$(BUILDDIR)/%_cu.o : %.cc *.h ../../src/*.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -Xcompiler -fPIC -c -x cu $< -o $@
endif

# Generic target and build rules: objects from C++ compilation
# (NB do not include CUINC here! add it only for NVTX or curand #679)
$(BUILDDIR)/%.o : %.cc *.h ../../src/*.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c $< -o $@

# Apply special build flags only to CrossSectionKernel.cc and gCrossSectionKernel.cu (no fast math, see #117 and #516)
ifeq ($(shell $(CXX) --version | grep ^nvc++),)
$(BUILDDIR)/CrossSectionKernels.o: CXXFLAGS := $(filter-out -ffast-math,$(CXXFLAGS))
$(BUILDDIR)/CrossSectionKernels.o: CXXFLAGS += -fno-fast-math
ifneq ($(NVCC),)
$(BUILDDIR)/gCrossSectionKernels.o: CUFLAGS += -Xcompiler -fno-fast-math
endif
endif

# Apply special build flags only to check_sa.o and gcheck_sa.o (NVTX in timermap.h, #679)
$(BUILDDIR)/check_sa.o: CXXFLAGS += $(USE_NVTX) $(CUINC)
$(BUILDDIR)/gcheck_sa.o: CXXFLAGS += $(USE_NVTX) $(CUINC)

# Apply special build flags only to check_sa and CurandRandomNumberKernel (curand headers, #679)
$(BUILDDIR)/check_sa.o: CXXFLAGS += $(CXXFLAGSCURAND)
$(BUILDDIR)/gcheck_sa.o: CUFLAGS += $(CXXFLAGSCURAND)
$(BUILDDIR)/CurandRandomNumberKernel.o: CXXFLAGS += $(CXXFLAGSCURAND)
$(BUILDDIR)/gCurandRandomNumberKernel.o: CUFLAGS += $(CXXFLAGSCURAND)
ifeq ($(RNDGEN),hasCurand)
$(BUILDDIR)/CurandRandomNumberKernel.o: CXXFLAGS += $(CUINC)
endif

# Avoid "warning: builtin __has_trivial_... is deprecated; use __is_trivially_... instead" in nvcc with icx2023 (#592)
ifneq ($(shell $(CXX) --version | egrep '^(Intel)'),)
ifneq ($(NVCC),)
CUFLAGS += -Xcompiler -Wno-deprecated-builtins
endif
endif

# Avoid clang warning "overriding '-ffp-contract=fast' option with '-ffp-contract=on'" (#516)
# This patch does remove the warning, but I prefer to keep it disabled for the moment...
###ifneq ($(shell $(CXX) --version | egrep '^(clang|Apple clang|Intel)'),)
###$(BUILDDIR)/CrossSectionKernels.o: CXXFLAGS += -Wno-overriding-t-option
###ifneq ($(NVCC),)
###$(BUILDDIR)/gCrossSectionKernels.o: CUFLAGS += -Xcompiler -Wno-overriding-t-option
###endif
###endif

#### Apply special build flags only to CPPProcess.cc (-flto)
###$(BUILDDIR)/CPPProcess.o: CXXFLAGS += -flto

#### Apply special build flags only to CPPProcess.cc (AVXFLAGS)
###$(BUILDDIR)/CPPProcess.o: CXXFLAGS += $(AVXFLAGS)

#-------------------------------------------------------------------------------

# Target (and build rules): common (src) library
commonlib : $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc $(BUILDDIR)/.build.$(TAG)
	$(MAKE) -C ../../src $(MAKEDEBUG) -f $(CUDACPP_SRC_MAKEFILE)

#-------------------------------------------------------------------------------

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')
###$(info processid_short=$(processid_short))

MG5AMC_CXXLIB = mg5amc_$(processid_short)_cpp
cxx_objects_lib=$(BUILDDIR)/CPPProcess.o $(BUILDDIR)/MatrixElementKernels.o $(BUILDDIR)/BridgeKernels.o $(BUILDDIR)/CrossSectionKernels.o
cxx_objects_exe=$(BUILDDIR)/CommonRandomNumberKernel.o $(BUILDDIR)/RamboSamplingKernels.o

ifneq ($(NVCC),)
MG5AMC_CULIB = mg5amc_$(processid_short)_cuda
cu_objects_lib=$(BUILDDIR)/gCPPProcess.o $(BUILDDIR)/gMatrixElementKernels.o $(BUILDDIR)/gBridgeKernels.o $(BUILDDIR)/gCrossSectionKernels.o
cu_objects_exe=$(BUILDDIR)/gCommonRandomNumberKernel.o $(BUILDDIR)/gRamboSamplingKernels.o
endif

# Target (and build rules): C++ and CUDA shared libraries
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(BUILDDIR)/fbridge.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: cxx_objects_lib += $(BUILDDIR)/fbridge.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib)
	$(CXX) -shared -o $@ $(cxx_objects_lib) $(CXXLIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

ifneq ($(NVCC),)
$(LIBDIR)/lib$(MG5AMC_CULIB).so: $(BUILDDIR)/fbridge_cu.o
$(LIBDIR)/lib$(MG5AMC_CULIB).so: cu_objects_lib += $(BUILDDIR)/fbridge_cu.o
$(LIBDIR)/lib$(MG5AMC_CULIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_objects_lib)
	$(NVCC) --shared -o $@ $(cu_objects_lib) $(CULIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)
endif

#-------------------------------------------------------------------------------

# Target (and build rules): Fortran include files
###$(INCDIR)/%.inc : ../%.inc
###	@if [ ! -d $(INCDIR) ]; then echo "mkdir -p $(INCDIR)"; mkdir -p $(INCDIR); fi
###	\cp $< $@

#-------------------------------------------------------------------------------

# Target (and build rules): C++ and CUDA standalone executables
$(cxx_main): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cxx_main): $(BUILDDIR)/check_sa.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel.o
	$(CXX) -o $@ $(BUILDDIR)/check_sa.o $(OMPFLAGS) -ldl -pthread $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CXXLIB) $(cxx_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel.o $(CURANDLIBFLAGS)

ifneq ($(NVCC),)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(cu_main): LIBFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(cu_main): LIBFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(cu_main): LIBFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif
$(cu_main): LIBFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cu_main): $(BUILDDIR)/gcheck_sa.o $(LIBDIR)/lib$(MG5AMC_CULIB).so $(cu_objects_exe) $(BUILDDIR)/gCurandRandomNumberKernel.o
	$(NVCC) -o $@ $(BUILDDIR)/gcheck_sa.o $(CUARCHFLAGS) $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CULIB) $(cu_objects_exe) $(BUILDDIR)/gCurandRandomNumberKernel.o $(CURANDLIBFLAGS)
endif

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

ifeq ($(UNAME_S),Darwin)
$(fcxx_main): LIBFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(fcxx_main): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(fcxx_main): $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe)
	$(CXX) -o $@ $(BUILDDIR)/fcheck_sa.o $(OMPFLAGS) $(BUILDDIR)/fsampler.o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_CXXLIB) $(cxx_objects_exe)

ifneq ($(NVCC),)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(fcu_main): LIBFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(fcu_main): LIBFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
endif
ifeq ($(UNAME_S),Darwin)
$(fcu_main): LIBFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(fcu_main): LIBFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(fcu_main): $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler_cu.o $(LIBDIR)/lib$(MG5AMC_CULIB).so $(cu_objects_exe)
	$(NVCC) -o $@ $(BUILDDIR)/fcheck_sa.o $(BUILDDIR)/fsampler_cu.o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_CULIB) $(cu_objects_exe)
endif

#-------------------------------------------------------------------------------

# Target (and build rules): test objects and test executable
$(BUILDDIR)/testxxx.o: $(GTESTLIBS)
$(BUILDDIR)/testxxx.o: INCFLAGS += $(GTESTINC)
$(BUILDDIR)/testxxx.o: testxxx_cc_ref.txt
$(testmain): $(BUILDDIR)/testxxx.o
$(testmain): cxx_objects_exe += $(BUILDDIR)/testxxx.o # Comment out this line to skip the C++ test of xxx functions

ifneq ($(NVCC),)
$(BUILDDIR)/testxxx_cu.o: $(GTESTLIBS)
$(BUILDDIR)/testxxx_cu.o: INCFLAGS += $(GTESTINC)
$(BUILDDIR)/testxxx_cu.o: testxxx_cc_ref.txt
$(testmain): $(BUILDDIR)/testxxx_cu.o
$(testmain): cu_objects_exe += $(BUILDDIR)/testxxx_cu.o # Comment out this line to skip the CUDA test of xxx functions
endif

$(BUILDDIR)/testmisc.o: $(GTESTLIBS)
$(BUILDDIR)/testmisc.o: INCFLAGS += $(GTESTINC)
$(testmain): $(BUILDDIR)/testmisc.o
$(testmain): cxx_objects_exe += $(BUILDDIR)/testmisc.o # Comment out this line to skip the C++ miscellaneous tests

ifneq ($(NVCC),)
$(BUILDDIR)/testmisc_cu.o: $(GTESTLIBS)
$(BUILDDIR)/testmisc_cu.o: INCFLAGS += $(GTESTINC)
$(testmain): $(BUILDDIR)/testmisc_cu.o
$(testmain): cu_objects_exe += $(BUILDDIR)/testmisc_cu.o # Comment out this line to skip the CUDA miscellaneous tests
endif

$(BUILDDIR)/runTest.o: $(GTESTLIBS)
$(BUILDDIR)/runTest.o: INCFLAGS += $(GTESTINC)
$(testmain): $(BUILDDIR)/runTest.o
$(testmain): cxx_objects_exe += $(BUILDDIR)/runTest.o

ifneq ($(NVCC),)
$(BUILDDIR)/runTest_cu.o: $(GTESTLIBS)
$(BUILDDIR)/runTest_cu.o: INCFLAGS += $(GTESTINC)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(testmain): LIBFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
$(testmain): LIBFLAGS += -lsvml # compile with icpx and link with nvcc (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(testmain): LIBFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif
$(testmain): $(BUILDDIR)/runTest_cu.o
$(testmain): cu_objects_exe  += $(BUILDDIR)/runTest_cu.o
endif

$(testmain): $(GTESTLIBS)
$(testmain): INCFLAGS +=  $(GTESTINC)
$(testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main

ifneq ($(OMPFLAGS),)
ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
$(testmain): LIBFLAGS += -liomp5 # see #578 (not '-qopenmp -static-intel' as in https://stackoverflow.com/questions/45909648)
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
$(testmain): LIBFLAGS += -L $(shell dirname $(shell $(CXX) -print-file-name=libc++.so)) -lomp # see #604
###else ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
###$(testmain): LIBFLAGS += ???? # OMP is not supported yet by cudacpp for Apple clang (see #578 and #604)
else
$(testmain): LIBFLAGS += -lgomp
endif
endif

ifeq ($(NVCC),) # link only runTest.o
$(testmain): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib) $(cxx_objects_exe) $(GTESTLIBS)
	$(CXX) -o $@ $(cxx_objects_lib) $(cxx_objects_exe) -ldl -pthread $(LIBFLAGS)
else # link both runTest.o and runTest_cu.o
$(testmain): LIBFLAGS += $(CULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib) $(cxx_objects_exe) $(cu_objects_lib) $(cu_objects_exe) $(GTESTLIBS)
	$(NVCC) -o $@ $(cxx_objects_lib) $(cxx_objects_exe) $(cu_objects_lib) $(cu_objects_exe) -ldl $(LIBFLAGS) -lcuda
endif

# Use target gtestlibs to build only googletest
ifneq ($(GTESTLIBS),)
gtestlibs: $(GTESTLIBS)
endif

# Use flock (Linux only, no Mac) to allow 'make -j' if googletest has not yet been downloaded https://stackoverflow.com/a/32666215
$(GTESTLIBS):
ifneq ($(shell which flock 2>/dev/null),)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	flock $(BUILDDIR)/.make_test.lock $(MAKE) -C $(TESTDIR)
else
	if [ -d $(TESTDIR) ]; then $(MAKE) -C $(TESTDIR); fi
endif

#-------------------------------------------------------------------------------

# Target: build all targets in all AVX modes (each AVX mode in a separate build directory)
# Split the avxall target into five separate targets to allow parallel 'make -j avxall' builds
# (Hack: add a fbridge.inc dependency to avxall, to ensure it is only copied once for all AVX modes)
avxnone:
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=none -f $(CUDACPP_MAKEFILE)

avxsse4:
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=sse4 -f $(CUDACPP_MAKEFILE)

avxavx2:
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=avx2 -f $(CUDACPP_MAKEFILE)

avx512y:
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512y -f $(CUDACPP_MAKEFILE)

avx512z:
	@echo
	$(MAKE) USEBUILDDIR=1 AVX=512z -f $(CUDACPP_MAKEFILE)

ifeq ($(UNAME_P),ppc64le)
###avxall: $(INCDIR)/fbridge.inc avxnone avxsse4
avxall: avxnone avxsse4
else ifeq ($(UNAME_P),arm)
###avxall: $(INCDIR)/fbridge.inc avxnone avxsse4
avxall: avxnone avxsse4
else
###avxall: $(INCDIR)/fbridge.inc avxnone avxsse4 avxavx2 avx512y avx512z
avxall: avxnone avxsse4 avxavx2 avx512y avx512z
endif

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
ifeq ($(USEBUILDDIR),1)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
	rm -f $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(LIBDIR)/lib$(MG5AMC_CULIB).so
endif
	$(MAKE) -C ../../src clean -f $(CUDACPP_SRC_MAKEFILE)
###	rm -rf $(INCDIR)

cleanall:
	@echo
	$(MAKE) USEBUILDDIR=0 clean -f $(CUDACPP_MAKEFILE)
	@echo
	$(MAKE) USEBUILDDIR=0 -C ../../src cleanall -f $(CUDACPP_SRC_MAKEFILE)
	rm -rf build.*

# Target: clean the builds as well as the gtest installation(s)
distclean: cleanall
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

# Target: check (run the C++ test executable)
# [NB THIS IS WHAT IS USED IN THE GITHUB CI!]
ifneq ($(NVCC),)
check: runTest cmpFcheck cmpFGcheck
else
check: runTest cmpFcheck
endif

# Target: runTest (run the C++ test executable runTest.exe)
runTest: all.$(TAG)
	$(RUNTIME) $(BUILDDIR)/runTest.exe

# Target: runCheck (run the C++ standalone executable check.exe, with a small number of events)
runCheck: all.$(TAG)
	$(RUNTIME) $(BUILDDIR)/check.exe -p 2 32 2

# Target: runGcheck (run the CUDA standalone executable gcheck.exe, with a small number of events)
runGcheck: all.$(TAG)
	$(RUNTIME) $(BUILDDIR)/gcheck.exe -p 2 32 2

# Target: runFcheck (run the Fortran standalone executable - with C++ MEs - fcheck.exe, with a small number of events)
runFcheck: all.$(TAG)
	$(RUNTIME) $(BUILDDIR)/fcheck.exe 2 32 2

# Target: runFGcheck (run the Fortran standalone executable - with CUDA MEs - fgcheck.exe, with a small number of events)
runFGcheck: all.$(TAG)
	$(RUNTIME) $(BUILDDIR)/fgcheck.exe 2 32 2

# Target: cmpFcheck (compare ME results from the C++ and Fortran with C++ MEs standalone executables, with a small number of events)
cmpFcheck: all.$(TAG)
	@echo
	@echo "$(BUILDDIR)/check.exe --common -p 2 32 2"
	@echo "$(BUILDDIR)/fcheck.exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(BUILDDIR)/check.exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(BUILDDIR)/fcheck.exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/C++)    = $${me1}"; echo "Avg ME (F77/C++)    = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/C++) returned NaN"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/C++) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi

# Target: cmpFGcheck (compare ME results from the CUDA and Fortran with CUDA MEs standalone executables, with a small number of events)
cmpFGcheck: all.$(TAG)
	@echo
	@echo "$(BUILDDIR)/gcheck.exe --common -p 2 32 2"
	@echo "$(BUILDDIR)/fgcheck.exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(BUILDDIR)/gcheck.exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(BUILDDIR)/fgcheck.exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/CUDA)   = $${me1}"; echo "Avg ME (F77/CUDA)   = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/CUDA) crashed"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/CUDA) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi

# Target: memcheck (run the CUDA standalone executable gcheck.exe with a small number of events through cuda-memcheck)
memcheck: all.$(TAG)
	$(RUNTIME) $(CUDA_HOME)/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes $(BUILDDIR)/gcheck.exe -p 2 32 2

#-------------------------------------------------------------------------------
