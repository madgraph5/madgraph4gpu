# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: S. Hageboeck, D. Massaro, O. Mattelaer, S. Roiser, J. Teig, A. Valassi (2020-2025) for the MG5aMC CUDACPP plugin.

#=== Determine the name of this makefile (https://ftp.gnu.org/old-gnu/Manuals/make-3.80/html_node/make_17.html)
#=== NB: use ':=' to ensure that the value of CUDACPP_MAKEFILE is not modified further down after including make_opts
#=== NB: use 'override' to ensure that the value can not be modified from the outside
override CUDACPP_MAKEFILE := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
###$(info CUDACPP_MAKEFILE='$(CUDACPP_MAKEFILE)')

#=== NB: different names (e.g. cudacpp.mk and cudacpp_src.mk) are used in the Subprocess and src directories
override CUDACPP_SRC_MAKEFILE = cudacpp_src.mk

#-------------------------------------------------------------------------------

#=== Include cudacpp_config.mk

# Check that the user-defined choices of BACKEND, FPTYPE, HELINL, HRDCOD are supported (and configure defaults if no user-defined choices exist)
# Stop with an error if BACKEND=cuda and nvcc is missing or if BACKEND=hip and hipcc is missing.
# Determine CUDACPP_BUILDDIR from a DIRTAG based on BACKEND, FPTYPE, HELINL, HRDCOD and from the user-defined choice of USEBUILDDIR
include ../../src/cudacpp_config.mk

# Export CUDACPP_BUILDDIR (so that there is no need to check/define it again in cudacpp_src.mk)
export CUDACPP_BUILDDIR

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

# OM: including make_opts is crucial for MG5aMC flag consistency/documentation
# AV: disable the inclusion of make_opts if the file has not been generated (standalone cudacpp)
ifneq ($(wildcard ../../Source/make_opts),)
  include ../../Source/make_opts
endif

#-------------------------------------------------------------------------------

#=== Redefine BACKEND if the current value is 'cppauto'

# Set the default BACKEND choice corresponding to 'cppauto' (the 'best' C++ vectorization available: eventually use native instead?)
ifeq ($(BACKEND),cppauto)
  ifeq ($(UNAME_P),ppc64le)
    override BACKEND = cppsse4
  else ifeq ($(UNAME_P),arm)
    override BACKEND = cppsse4
  else ifeq ($(wildcard /proc/cpuinfo),)
    override BACKEND = cppnone
    ###$(warning Using BACKEND='$(BACKEND)' because host SIMD features cannot be read from /proc/cpuinfo)
  else ifeq ($(shell grep -m1 -c avx512vl /proc/cpuinfo)$(shell $(CXX) --version | grep ^clang),1)
    override BACKEND = cpp512y
  else
    override BACKEND = cppavx2
    ###ifneq ($(shell grep -m1 -c avx512vl /proc/cpuinfo),1)
    ###  $(warning Using BACKEND='$(BACKEND)' because host does not support avx512vl)
    ###else
    ###  $(warning Using BACKEND='$(BACKEND)' because this is faster than avx512vl for clang)
    ###endif
  endif
  $(info BACKEND=$(BACKEND) (was cppauto))
else
  $(info BACKEND='$(BACKEND)')
endif

#-------------------------------------------------------------------------------

#=== Configure the C++ compiler

CXXFLAGS = $(OPTFLAGS) -std=c++17 -Wall -Wshadow -Wextra
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

# Export CXXFLAGS (so that there is no need to check/define it again in cudacpp_src.mk)
export CXXFLAGS

#-------------------------------------------------------------------------------

#=== Configure the GPU compiler (CUDA or HIP)
#=== (note, this is done also for C++, as NVTX and CURAND/ROCRAND are also needed by the C++ backends)

# Set CUDA_HOME from the path to nvcc, if it exists
override CUDA_HOME = $(patsubst %/bin/nvcc,%,$(shell which nvcc 2>/dev/null))

# Set HIP_HOME from the path to hipcc, if it exists
override HIP_HOME = $(shell hipconfig --rocmpath)

# Configure CUDA_INC (for CURAND and NVTX) and NVTX if a CUDA installation exists (see #965)
ifeq ($(CUDA_HOME),)
  # CUDA_HOME is empty (nvcc not found)
  override CUDA_INC=
else ifeq ($(wildcard $(CUDA_HOME)/include/),)
  # CUDA_HOME is defined (nvcc was found) but $(CUDA_HOME)/include/ does not exist?
  override CUDA_INC=
else
  CUDA_INC = -I$(CUDA_HOME)/include/
endif
###$(info CUDA_INC=$(CUDA_INC))

# Configure NVTX if a CUDA include directory exists and NVTX headers exist (see #965)
ifeq ($(CUDA_INC),)
  # $(CUDA_HOME)/include/ does not exist
  override USE_NVTX=
else ifeq ($(wildcard $(CUDA_HOME)/include/nvtx3/nvToolsExt.h),)
  # $(CUDA_HOME)/include/ exists but NVTX headers do not exist?
  override USE_NVTX=
else
  # $(CUDA_HOME)/include/nvtx.h exists: use NVTX
  # (NB: the option to disable NVTX if 'USE_NVTX=' is defined has been removed)
  override USE_NVTX=-DUSE_NVTX
endif
###$(info USE_NVTX=$(USE_NVTX))

# NB: NEW LOGIC FOR ENABLING AND DISABLING CUDA OR HIP BUILDS (AV Feb-Mar 2024)
# - In the old implementation, by default the C++ targets for one specific AVX were always built together with either CUDA or HIP.
# If both CUDA and HIP were installed, then CUDA took precedence over HIP, and the only way to force HIP builds was to disable
# CUDA builds by setting CUDA_HOME to an invalid value (as CUDA_HOME took precdence over PATH to find the installation of nvcc).
# Similarly, C++-only builds could be forced by setting CUDA_HOME and/or HIP_HOME to invalid values. A check for an invalid nvcc
# in CUDA_HOME or an invalid hipcc HIP_HOME was necessary to ensure this logic, and had to be performed at the very beginning.
# - In the new implementation (PR #798), separate individual builds are performed for one specific C++/AVX mode, for CUDA or
# for HIP. The choice of the type of build is taken depending on the value of the BACKEND variable (replacing the AVX variable).
# Unlike what happened in the past, nvcc and hipcc must have already been added to PATH. Using 'which nvcc' and 'which hipcc',
# their existence and their location is checked, and the variables CUDA_HOME and HIP_HOME are internally set by this makefile.
# This must be still done before backend-specific customizations, e.g. because CURAND and NVTX are also used in C++ builds.
# Note also that a preliminary check for nvcc and hipcc if BACKEND is cuda or hip is performed in cudacpp_config.mk.
# - Note also that the REQUIRE_CUDA variable (which was used in the past, e.g. for CI tests on GPU #443) is now (PR #798) no
# longer necessary, as it is now equivalent to BACKEND=cuda. Similarly, there is no need to introduce a REQUIRE_HIP variable.

#=== Configure the CUDA or HIP compiler (only for the CUDA and HIP backends)
#=== (NB: throughout all makefiles, an empty GPUCC is used to indicate that this is a C++ build, i.e. that BACKEND is neither cuda nor hip!)

ifeq ($(BACKEND),cuda)

  # If CXX is not a single word (example "clang++ --gcc-toolchain...") then disable CUDA builds (issue #505)
  # This is because it is impossible to pass this to "GPUFLAGS += -ccbin <host-compiler>" below
  ifneq ($(words $(subst ccache ,,$(CXX))),1) # allow at most "CXX=ccache <host-compiler>" from outside
    $(error BACKEND=$(BACKEND) but CUDA builds are not supported for multi-word CXX "$(CXX)")
  endif

  # Set GPUCC as $(CUDA_HOME)/bin/nvcc (it was already checked above that this exists)
  GPUCC = $(CUDA_HOME)/bin/nvcc
  XCOMPILERFLAG = -Xcompiler
  GPULANGUAGE = cu
  GPUSUFFIX = cuda

  # Optimization flags
  GPUFLAGS = $(foreach opt, $(OPTFLAGS), $(XCOMPILERFLAG) $(opt))

  # NVidia CUDA architecture flags
  # See https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html
  # See https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
  # Default: use compute capability 70 for V100 (CERN lxbatch, CERN itscrd, Juwels Cluster).
  # This will embed device code for 70, and PTX for 70+.
  # One may pass MADGRAPH_CUDA_ARCHITECTURE (comma-separated list) to the make command to use another value or list of values (see #533).
  # Examples: use 60 for P100 (Piz Daint), 80 for A100 (Juwels Booster, NVidia raplab/Curiosity).
  MADGRAPH_CUDA_ARCHITECTURE ?= 70
  ###GPUARCHFLAGS = -gencode arch=compute_$(MADGRAPH_CUDA_ARCHITECTURE),code=compute_$(MADGRAPH_CUDA_ARCHITECTURE) -gencode arch=compute_$(MADGRAPH_CUDA_ARCHITECTURE),code=sm_$(MADGRAPH_CUDA_ARCHITECTURE) # Older implementation (AV): go back to this one for multi-GPU support #533
  ###GPUARCHFLAGS = --gpu-architecture=compute_$(MADGRAPH_CUDA_ARCHITECTURE) --gpu-code=sm_$(MADGRAPH_CUDA_ARCHITECTURE),compute_$(MADGRAPH_CUDA_ARCHITECTURE)  # Newer implementation (SH): cannot use this as-is for multi-GPU support #533
  comma:=,
  GPUARCHFLAGS = $(foreach arch,$(subst $(comma), ,$(MADGRAPH_CUDA_ARCHITECTURE)),-gencode arch=compute_$(arch),code=compute_$(arch) -gencode arch=compute_$(arch),code=sm_$(arch))
  GPUFLAGS += $(GPUARCHFLAGS)

  # Other NVidia-specific flags
  CUDA_OPTFLAGS = -lineinfo
  GPUFLAGS += $(CUDA_OPTFLAGS)

  # NVCC version
  ###GPUCC_VERSION = $(shell $(GPUCC) --version | grep 'Cuda compilation tools' | cut -d' ' -f5 | cut -d, -f1)

  # Fast math
  GPUFLAGS += -use_fast_math

  # Extra build warnings
  GPUFLAGS += $(XCOMPILERFLAG) -Wunused-parameter
  ###GPUFLAGS += $(XCOMPILERFLAG) -Wall $(XCOMPILERFLAG) -Wextra $(XCOMPILERFLAG) -Wshadow

  # CUDA includes and NVTX
  GPUFLAGS += $(CUDA_INC) $(USE_NVTX) 

  # C++ standard
  GPUFLAGS += -std=c++17 # need CUDA >= 11.2 (see #333): this is enforced in mgOnGpuConfig.h

  # For nvcc, use -maxrregcount to control the maximum number of registries (this does not exist in hipcc)
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###GPUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###GPUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###GPUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###GPUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)

  # Set the host C++ compiler for nvcc via "-ccbin <host-compiler>"
  # (NB issue #505: this must be a single word, "clang++ --gcc-toolchain..." is not supported)
  GPUFLAGS += -ccbin $(shell which $(subst ccache ,,$(CXX)))

  # Allow newer (unsupported) C++ compilers with older versions of CUDA if ALLOW_UNSUPPORTED_COMPILER_IN_CUDA is set (#504)
  ifneq ($(origin ALLOW_UNSUPPORTED_COMPILER_IN_CUDA),undefined)
    GPUFLAGS += -allow-unsupported-compiler
  endif

else ifeq ($(BACKEND),hip)

  # example architecture values MI200:gfx90a, MI350X:gfx942
  MADGRAPH_HIP_ARCHITECTURE ?= gfx942
  # Set GPUCC as $(HIP_HOME)/bin/hipcc (it was already checked above that this exists)
  GPUCC = $(HIP_HOME)/bin/hipcc
  XCOMPILERFLAG =
  GPULANGUAGE = hip
  GPUSUFFIX = hip

  # Optimization flags
  override OPTFLAGS = -O2 # work around "Memory access fault" in gq_ttq for HIP #806: disable hipcc -O3 optimizations
  GPUFLAGS = $(foreach opt, $(OPTFLAGS), $(XCOMPILERFLAG) $(opt))

  # DEBUG FLAGS (for #806: see https://hackmd.io/@gmarkoma/lumi_finland)
  ###GPUFLAGS += -ggdb # FOR DEBUGGING ONLY

  # AMD HIP architecture flags
  GPUARCHFLAGS = --offload-arch=${MADGRAPH_HIP_ARCHITECTURE}
  GPUFLAGS += $(GPUARCHFLAGS)

  # Other AMD-specific flags
  GPUFLAGS += -target x86_64-linux-gnu -DHIP_PLATFORM=amd

  # Fast math (is -DHIP_FAST_MATH equivalent to -ffast-math?)
  GPUFLAGS += -DHIP_FAST_MATH

  # Extra build warnings
  ###GPUFLAGS += $(XCOMPILERFLAG) -Wall $(XCOMPILERFLAG) -Wextra $(XCOMPILERFLAG) -Wshadow

  # HIP includes
  HIP_INC = -I$(HIP_HOME)/include/
  GPUFLAGS += $(HIP_INC)

  # C++ standard
  GPUFLAGS += -std=c++17

else

  # Backend is neither cuda nor hip
  override GPUCC=
  override GPUFLAGS=

  # Sanity check, this should never happen: if GPUCC is empty, then this is a C++ build, i.e. BACKEND is neither cuda nor hip.
  # In practice, in the following, "ifeq ($(GPUCC),)" is equivalent to "ifneq ($(findstring cpp,$(BACKEND)),)".
  # Conversely, note that GPUFLAGS is non-empty also for C++ builds, but it is never used in that case.
  ifeq ($(findstring cpp,$(BACKEND)),)
    $(error INTERNAL ERROR! Unknown backend BACKEND='$(BACKEND)': supported backends are $(foreach backend,$(SUPPORTED_BACKENDS),'$(backend)'))
  endif

endif

# Export GPUCC, GPUFLAGS, GPULANGUAGE, GPUSUFFIX (so that there is no need to check/define them again in cudacpp_src.mk)
export GPUCC
export GPUFLAGS
export GPULANGUAGE
export GPUSUFFIX

#-------------------------------------------------------------------------------

#=== Configure ccache for C++ and CUDA/HIP builds

# Enable ccache only if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif
#ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
#  override AR:=ccache $(AR)
#endif
ifneq ($(GPUCC),)
  ifeq ($(USECCACHE)$(shell echo $(GPUCC) | grep ccache),1)
    override GPUCC:=ccache $(GPUCC)
  endif
endif

#-------------------------------------------------------------------------------

#=== Configure common compiler flags for C++ and CUDA/HIP

INCFLAGS = -I.
OPTFLAGS = -O3 # this ends up in GPUFLAGS too (should it?), cannot add -Ofast or -ffast-math here

# Dependency on src directory
ifeq ($(GPUCC),)
MG5AMC_COMMONLIB = mg5amc_common_cpp
else
MG5AMC_COMMONLIB = mg5amc_common_$(GPUSUFFIX)
endif
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

# Export CXXNAMESUFFIX (so that there is no need to check/define it again in cudacpp_test.mk)
export CXXNAMESUFFIX

# Dependency on test directory
# Within the madgraph4gpu git repo: by default use a common gtest installation in <topdir>/test (optionally use an external or local gtest)
# Outside the madgraph4gpu git repo: by default do not build the tests (optionally use an external or local gtest)
# Do not build the tests if USEGTEST is equal to 0 (default inside launch_plugin.py, see https://github.com/madgraph5/madgraph4gpu/issues/878)
###GTEST_ROOT = /cvmfs/sft.cern.ch/lcg/releases/gtest/1.11.0-21e8c/x86_64-centos8-gcc11-opt/# example of an external gtest installation
###LOCALGTEST = yes# comment this out (or use make LOCALGTEST=yes) to build tests using a local gtest installation
TESTDIRCOMMON = ../../../../../test
TESTDIRLOCAL = ../../test
ifeq ($(USEGTEST),0)
  TESTDIR=
  GTEST_ROOT=
else ifneq ($(wildcard $(GTEST_ROOT)),)
  TESTDIR=
else ifneq ($(LOCALGTEST),)
  TESTDIR=$(TESTDIRLOCAL)
  GTEST_ROOT=$(TESTDIR)/googletest/install$(CXXNAMESUFFIX)
else ifneq ($(wildcard ../../../../../epochX/cudacpp/CODEGEN),)
  TESTDIR=$(TESTDIRCOMMON)
  GTEST_ROOT= $(TESTDIR)/googletest/install$(CXXNAMESUFFIX)
else
  TESTDIR=
endif
ifneq ($(GTEST_ROOT),)
  GTESTLIBDIR = $(GTEST_ROOT)/lib64/
  GTESTLIBS = $(GTESTLIBDIR)/libgtest.a
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

#=== Configure PowerPC-specific compiler flags for C++ and CUDA/HIP

# PowerPC-specific CXX compiler flags (being reviewed)
ifeq ($(UNAME_P),ppc64le)
  CXXFLAGS+= -mcpu=power9 -mtune=power9 # gains ~2-3% both for cppnone and cppsse4
  # Throughput references without the extra flags below: cppnone=1.41-1.42E6, cppsse4=2.15-2.19E6
  ###CXXFLAGS+= -DNO_WARN_X86_INTRINSICS # no change
  ###CXXFLAGS+= -fpeel-loops # no change
  ###CXXFLAGS+= -funroll-loops # gains ~1% for cppnone, loses ~1% for cppsse4
  ###CXXFLAGS+= -ftree-vectorize # no change
  ###CXXFLAGS+= -flto # would increase to cppnone=4.08-4.12E6, cppsse4=4.99-5.03E6!
else
  ###CXXFLAGS+= -flto # also on Intel this would increase throughputs by a factor 2 to 4...
  ######CXXFLAGS+= -fno-semantic-interposition # no benefit (neither alone, nor combined with -flto)
endif

# PowerPC-specific CUDA/HIP compiler flags (to be reviewed!)
ifeq ($(UNAME_P),ppc64le)
  GPUFLAGS+= $(XCOMPILERFLAG) -mno-float128
endif

#-------------------------------------------------------------------------------

#=== Configure defaults for OMPFLAGS

# Disable OpenMP by default: enable OpenMP only if USEOPENMP=1 (#758)
ifeq ($(USEOPENMP),1)
  ###$(info USEOPENMP==1: will build with OpenMP if possible)
  ifneq ($(findstring hipcc,$(GPUCC)),)
    override OMPFLAGS = # disable OpenMP MT when using hipcc #802
  else ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
    override OMPFLAGS = -fopenmp
    ###override OMPFLAGS = # disable OpenMP MT on Intel (was ok without GPUCC but not ok with GPUCC before #578)
  else ifneq ($(shell $(CXX) --version | egrep '^clang version 16'),)
    ###override OMPFLAGS = # disable OpenMP on clang16 #904
    $(error OpenMP is not supported by cudacpp on clang16 - issue #904)
  else ifneq ($(shell $(CXX) --version | egrep '^clang version 17'),)
    ###override OMPFLAGS = # disable OpenMP on clang17 #904
    $(error OpenMP is not supported by cudacpp on clang17 - issue #904)
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
else
  ###$(info USEOPENMP!=1: will build without OpenMP)
  override OMPFLAGS =
endif

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for RNDGEN (legacy!), HASCURAND, HASHIPRAND

# If the legacy RNDGEN exists, this take precedence over any HASCURAND choice (but a warning is printed out)
###$(info RNDGEN=$(RNDGEN))
ifneq ($(RNDGEN),)
  $(warning Environment variable RNDGEN is no longer supported, please use HASCURAND instead!)
  ifeq ($(RNDGEN),hasCurand)
    override HASCURAND = $(RNDGEN)
  else ifeq ($(RNDGEN),hasNoCurand)
    override HASCURAND = $(RNDGEN)
  else ifneq ($(RNDGEN),hasNoCurand)
    $(error Unknown RNDGEN='$(RNDGEN)': only 'hasCurand' and 'hasNoCurand' are supported - but use HASCURAND instead!)
  endif
endif

# Set the default HASCURAND (curand random number generator) choice, if no prior choice exists for HASCURAND
# (NB: allow HASCURAND=hasCurand even if $(GPUCC) does not point to nvcc: assume CUDA_HOME was defined correctly...)
ifeq ($(HASCURAND),)
  ifeq ($(GPUCC),) # CPU-only build
    ifeq ($(CUDA_INC),)
      # $(CUDA_HOME)/include/ does not exist (see #965)
      override HASCURAND = hasNoCurand
    else ifeq ($(wildcard $(CUDA_HOME)/include/curand.h),)
      # $(CUDA_HOME)/include/ exists but CURAND headers do not exist? (see #965)
      override HASCURAND = hasNoCurand
    else
      # By default, assume that curand is installed if a CUDA installation exists
      override HASCURAND = hasCurand
    endif
  else ifeq ($(findstring nvcc,$(GPUCC)),nvcc) # Nvidia GPU build
    # By default, assume that curand is installed if a CUDA build is requested
    override HASCURAND = hasCurand
  else # non-Nvidia GPU build
    override HASCURAND = hasNoCurand
  endif
endif

# Set the default HASHIPRAND (hiprand random number generator) choice, if no prior choice exists for HASHIPRAND
# (NB: allow HASHIPRAND=hasHiprand even if $(GPUCC) does not point to hipcc: assume HIP_HOME was defined correctly...)
ifeq ($(HASHIPRAND),)
  ifeq ($(GPUCC),) # CPU-only build
    override HASHIPRAND = hasNoHiprand
  else ifeq ($(findstring hipcc,$(GPUCC)),hipcc) # AMD GPU build
    override HASHIPRAND = hasHiprand
  else # non-AMD GPU build
    override HASHIPRAND = hasNoHiprand
  endif
endif

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for HASBLAS

# Set the default HASBLAS (cuBLAS/hipBLAS) choice and check prior choices for HASBLAS

ifeq ($(HASBLAS),)
  ifeq ($(GPUCC),) # CPU-only build
    override HASBLAS = hasNoBlas
  else ifeq ($(findstring nvcc,$(GPUCC)),nvcc) # Nvidia GPU build
    ifeq ($(wildcard $(CUDA_HOME)/include/cublas_v2.h),)
      # cuBLAS headers do not exist??
      override HASBLAS = hasNoBlas
    else
      override HASBLAS = hasBlas
    endif
  else ifeq ($(findstring hipcc,$(GPUCC)),hipcc) # AMD GPU build
    ifeq ($(wildcard $(HIP_HOME)/include/hipblas/hipblas.h),)
      # hipBLAS headers do not exist??
      override HASBLAS = hasNoBlas
    else
      override HASBLAS = hasBlas
    endif
  else
    override HASBLAS = hasNoBlas
  endif
endif

#-------------------------------------------------------------------------------

#=== Set the CUDA/HIP/C++ compiler flags appropriate to user-defined choices of AVX, FPTYPE, HELINL, HRDCOD

# Set the build flags appropriate to OMPFLAGS
$(info OMPFLAGS=$(OMPFLAGS))
CXXFLAGS += $(OMPFLAGS)

# Set the build flags appropriate to each BACKEND choice (example: "make BACKEND=cppnone")
# [NB MGONGPU_PVW512 is needed because "-mprefer-vector-width=256" is not exposed in a macro]
# [See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96476]
ifeq ($(UNAME_P),ppc64le)
  ifeq ($(BACKEND),cppsse4)
    override AVXFLAGS = -D__SSE4_2__ # Power9 VSX with 128 width (VSR registers)
  else ifeq ($(BACKEND),cppavx2)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on PowerPC for the moment)
  else ifeq ($(BACKEND),cpp512y)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on PowerPC for the moment)
  else ifeq ($(BACKEND),cpp512z)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on PowerPC for the moment)
  endif
else ifeq ($(UNAME_P),arm)
  ifeq ($(BACKEND),cppsse4)
    override AVXFLAGS = -D__SSE4_2__ # ARM NEON with 128 width (Q/quadword registers)
  else ifeq ($(BACKEND),cppavx2)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on ARM for the moment)
  else ifeq ($(BACKEND),cpp512y)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on ARM for the moment)
  else ifeq ($(BACKEND),cpp512z)
    $(error Invalid SIMD BACKEND='$(BACKEND)': only 'cppnone' and 'cppsse4' are supported on ARM for the moment)
  endif
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
  ifeq ($(BACKEND),cppnone)
    override AVXFLAGS = -mno-sse3 # no SIMD
  else ifeq ($(BACKEND),cppsse4)
    override AVXFLAGS = -mno-avx # SSE4.2 with 128 width (xmm registers)
  else ifeq ($(BACKEND),cppavx2)
    override AVXFLAGS = -march=haswell # AVX2 with 256 width (ymm registers) [DEFAULT for clang]
  else ifeq ($(BACKEND),cpp512y)
    override AVXFLAGS = -march=skylake -mprefer-vector-width=256 # AVX512 with 256 width (ymm registers) [DEFAULT for gcc]
  else ifeq ($(BACKEND),cpp512z)
    override AVXFLAGS = -march=skylake -DMGONGPU_PVW512 # AVX512 with 512 width (zmm registers)
  endif
else
  ifeq ($(BACKEND),cppnone)
    override AVXFLAGS = -march=x86-64 # no SIMD (see #588)
  else ifeq ($(BACKEND),cppsse4)
    override AVXFLAGS = -march=nehalem # SSE4.2 with 128 width (xmm registers)
  else ifeq ($(BACKEND),cppavx2)
    override AVXFLAGS = -march=haswell # AVX2 with 256 width (ymm registers) [DEFAULT for clang]
  else ifeq ($(BACKEND),cpp512y)
    override AVXFLAGS = -march=skylake-avx512 -mprefer-vector-width=256 # AVX512 with 256 width (ymm registers) [DEFAULT for gcc]
  else ifeq ($(BACKEND),cpp512z)
    override AVXFLAGS = -march=skylake-avx512 -DMGONGPU_PVW512 # AVX512 with 512 width (zmm registers)
  endif
endif
# For the moment, use AVXFLAGS everywhere (in C++ builds): eventually, use them only in encapsulated implementations?
ifeq ($(GPUCC),)
  CXXFLAGS+= $(AVXFLAGS)
endif

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
$(info FPTYPE='$(FPTYPE)')
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_DOUBLE
  GPUFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT -DMGONGPU_FPTYPE2_FLOAT
  GPUFLAGS += -DMGONGPU_FPTYPE_FLOAT -DMGONGPU_FPTYPE2_FLOAT
else ifeq ($(FPTYPE),m)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_FLOAT
  GPUFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'd', 'f' and 'm' are supported)
endif

# Set the build flags appropriate to each HELINL choice (example: "make HELINL=1")
$(info HELINL='$(HELINL)')
ifeq ($(HELINL),1)
  CXXFLAGS += -DMGONGPU_INLINE_HELAMPS
  GPUFLAGS += -DMGONGPU_INLINE_HELAMPS
else ifneq ($(HELINL),0)
  $(error Unknown HELINL='$(HELINL)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each HRDCOD choice (example: "make HRDCOD=1")
$(info HRDCOD='$(HRDCOD)')
ifeq ($(HRDCOD),1)
  CXXFLAGS += -DMGONGPU_HARDCODE_PARAM
  GPUFLAGS += -DMGONGPU_HARDCODE_PARAM
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

#=== Set the CUDA/HIP/C++ compiler and linker flags appropriate to user-defined choices of HASCURAND, HASHIPRAND

$(info HASCURAND=$(HASCURAND))
$(info HASHIPRAND=$(HASHIPRAND))
override RNDCXXFLAGS=
override RNDLIBFLAGS=

# Set the RNDCXXFLAGS and RNDLIBFLAGS build flags appropriate to each HASCURAND choice (example: "make HASCURAND=hasNoCurand")
ifeq ($(HASCURAND),hasNoCurand)
  override RNDCXXFLAGS += -DMGONGPU_HAS_NO_CURAND
else ifeq ($(HASCURAND),hasCurand)
  override RNDLIBFLAGS += -L$(CUDA_HOME)/lib64/ -lcurand # NB: -lcuda is not needed here!
else
  $(error Unknown HASCURAND='$(HASCURAND)': only 'hasCurand' and 'hasNoCurand' are supported)
endif

# Set the RNDCXXFLAGS and RNDLIBFLAGS build flags appropriate to each HASHIPRAND choice (example: "make HASHIPRAND=hasNoHiprand")
ifeq ($(HASHIPRAND),hasNoHiprand)
  override RNDCXXFLAGS += -DMGONGPU_HAS_NO_HIPRAND
else ifeq ($(HASHIPRAND),hasHiprand)
  override RNDLIBFLAGS += -L$(HIP_HOME)/lib/ -lhiprand
else ifneq ($(HASHIPRAND),hasHiprand)
  $(error Unknown HASHIPRAND='$(HASHIPRAND)': only 'hasHiprand' and 'hasNoHiprand' are supported)
endif

#$(info RNDCXXFLAGS=$(RNDCXXFLAGS))
#$(info RNDLIBFLAGS=$(RNDLIBFLAGS))

#=== Set the CUDA/HIP/C++ compiler and linker flags appropriate to user-defined choices of HASBLAS

$(info HASBLAS=$(HASBLAS))
override BLASCXXFLAGS=
override BLASLIBFLAGS=

# Set the RNDCXXFLAGS and RNDLIBFLAGS build flags appropriate to each HASBLAS choice (example: "make HASBLAS=hasNoBlas")
ifeq ($(HASBLAS),hasNoBlas)
  override BLASCXXFLAGS += -DMGONGPU_HAS_NO_BLAS
else ifeq ($(HASBLAS),hasBlas)
  ifeq ($(findstring nvcc,$(GPUCC)),nvcc) # Nvidia GPU build
    override BLASLIBFLAGS = -L$(CUDA_HOME)/lib64/ -lcublas
  else ifeq ($(findstring hipcc,$(GPUCC)),hipcc) # AMD GPU build
    override BLASLIBFLAGS = -L$(HIP_HOME)/lib/ -lhipblas
  endif
else
  $(error Unknown HASBLAS='$(HASBLAS)': only 'hasBlas' and 'hasNoBlas' are supported)
endif
CXXFLAGS += $(BLASCXXFLAGS)
GPUFLAGS += $(BLASCXXFLAGS)

#$(info BLASCXXFLAGS=$(BLASCXXFLAGS))
#$(info BLASLIBFLAGS=$(BLASLIBFLAGS))

#-------------------------------------------------------------------------------

#=== Configure Position-Independent Code
CXXFLAGS += -fPIC
GPUFLAGS += $(XCOMPILERFLAG) -fPIC

#-------------------------------------------------------------------------------

#=== Configure channelid debugging
ifneq ($(MG5AMC_CHANNELID_DEBUG),)
  CXXFLAGS += -DMGONGPU_CHANNELID_DEBUG
  GPUFLAGS += -DMGONGPU_CHANNELID_DEBUG
endif

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of builds with different random number generators)
override TAG = $(patsubst cpp%,%,$(BACKEND))_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)_$(HASCURAND)_$(HASHIPRAND)

# Export TAG (so that there is no need to check/define it again in cudacpp_src.mk)
export TAG

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
override BUILDDIR = $(CUDACPP_BUILDDIR)
ifeq ($(USEBUILDDIR),1)
  override LIBDIR = ../../lib/$(BUILDDIR)
  override LIBDIRRPATH = '$$ORIGIN/../$(LIBDIR)'
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR == 1))
else
  override LIBDIR = ../../lib
  override LIBDIRRPATH = '$$ORIGIN/$(LIBDIR)'
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR != 1))
endif
###override INCDIR = ../../include
###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

# On Linux, set rpath to LIBDIR to make it unnecessary to use LD_LIBRARY_PATH
# Use relative paths with respect to the executables or shared libraries ($ORIGIN on Linux)
# On Darwin, building libraries with absolute paths in LIBDIR makes this unnecessary
ifeq ($(UNAME_S),Darwin)
  override CXXLIBFLAGSRPATH =
  override GPULIBFLAGSRPATH =
  override CXXLIBFLAGSRPATH2 =
  override GPULIBFLAGSRPATH2 =
else
  # RPATH to gpu/cpp libs when linking executables
  override CXXLIBFLAGSRPATH = -Wl,-rpath=$(LIBDIRRPATH)
  override GPULIBFLAGSRPATH = -Xlinker -rpath=$(LIBDIRRPATH)
  # RPATH to common lib when linking gpu/cpp libs
  override CXXLIBFLAGSRPATH2 = -Wl,-rpath='$$ORIGIN'
  override GPULIBFLAGSRPATH2 = -Xlinker -rpath='$$ORIGIN'
endif

# Setting LD_LIBRARY_PATH or DYLD_LIBRARY_PATH in the RUNTIME is no longer necessary (neither on Linux nor on Mac)
override RUNTIME =

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================


ifeq ($(GPUCC),)
  cxx_checkmain=$(BUILDDIR)/check_cpp.exe
  cxx_fcheckmain=$(BUILDDIR)/fcheck_cpp.exe
  cxx_testmain=$(BUILDDIR)/runTest_cpp.exe
else
  gpu_checkmain=$(BUILDDIR)/check_$(GPUSUFFIX).exe
  gpu_fcheckmain=$(BUILDDIR)/fcheck_$(GPUSUFFIX).exe
  gpu_testmain=$(BUILDDIR)/runTest_$(GPUSUFFIX).exe
endif

# Explicitly define the default goal (this is not necessary as it is the first target, which is implicitly the default goal)
.DEFAULT_GOAL := all.$(TAG)

# First target (default goal)
ifeq ($(GPUCC),)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_checkmain) $(cxx_fcheckmain) $(if $(GTESTLIBS),$(cxx_testmain))
else
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(gpu_checkmain) $(gpu_fcheckmain) $(if $(GTESTLIBS),$(gpu_testmain))
endif

# Target (and build options): debug
MAKEDEBUG=
debug: OPTFLAGS = -g -O0
debug: CUDA_OPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

# Target (and build options): address sanitizer #207
###CXXLIBFLAGSASAN =
###GPULIBFLAGSASAN =
###asan: OPTFLAGS = -g -O0 -fsanitize=address -fno-omit-frame-pointer
###asan: CUDA_OPTFLAGS = -G $(XCOMPILERFLAG) -fsanitize=address $(XCOMPILERFLAG) -fno-omit-frame-pointer
###asan: CXXLIBFLAGSASAN = -fsanitize=address
###asan: GPULIBFLAGSASAN = -Xlinker -fsanitize=address -Xlinker -shared
###asan: MAKEDEBUG := debug
###asan: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtagsb)" != "" ]; then echo "Cannot build for tag=$(TAG) as old builds exist for other tags:"; echo "  $(oldtagsb)"; echo "Please run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

# Apply special build flags only to CrossSectionKernel_<cpp|$(GPUSUFFIX)>.o (no fast math, see #117 and #516)
# Added edgecase for HIP compilation
ifeq ($(shell $(CXX) --version | grep ^nvc++),)
$(BUILDDIR)/CrossSectionKernels_cpp.o: CXXFLAGS := $(filter-out -ffast-math,$(CXXFLAGS))
$(BUILDDIR)/CrossSectionKernels_cpp.o: CXXFLAGS += -fno-fast-math
$(BUILDDIR)/CrossSectionKernels_$(GPUSUFFIX).o: GPUFLAGS += $(XCOMPILERFLAG) -fno-fast-math
endif

# Apply special build flags only to check_sa_<cpp|$(GPUSUFFIX)>.o (NVTX in timermap.h, #679)
$(BUILDDIR)/check_sa_cpp.o: CXXFLAGS += $(USE_NVTX) $(CUDA_INC)
$(BUILDDIR)/check_sa_$(GPUSUFFIX).o: CXXFLAGS += $(USE_NVTX) $(CUDA_INC)

# Apply special build flags only to check_sa_<cpp|$(GPUSUFFIX)>.o and (Cu|Hip)randRandomNumberKernel_<cpp|$(GPUSUFFIX)>.o
$(BUILDDIR)/check_sa_cpp.o: CXXFLAGS += $(RNDCXXFLAGS)
$(BUILDDIR)/check_sa_$(GPUSUFFIX).o: GPUFLAGS += $(RNDCXXFLAGS)
$(BUILDDIR)/CurandRandomNumberKernel_cpp.o: CXXFLAGS += $(RNDCXXFLAGS)
$(BUILDDIR)/CurandRandomNumberKernel_$(GPUSUFFIX).o: GPUFLAGS += $(RNDCXXFLAGS)
$(BUILDDIR)/HiprandRandomNumberKernel_cpp.o: CXXFLAGS += $(RNDCXXFLAGS)
$(BUILDDIR)/HiprandRandomNumberKernel_$(GPUSUFFIX).o: GPUFLAGS += $(RNDCXXFLAGS)
ifeq ($(HASCURAND),hasCurand) # curand headers, #679
$(BUILDDIR)/CurandRandomNumberKernel_cpp.o: CXXFLAGS += $(CUDA_INC)
endif
ifeq ($(HASHIPRAND),hasHiprand) # hiprand headers
$(BUILDDIR)/HiprandRandomNumberKernel_cpp.o: CXXFLAGS += $(HIP_INC)
endif

# Avoid "warning: builtin __has_trivial_... is deprecated; use __is_trivially_... instead" in GPUCC with icx2023 (#592)
ifneq ($(shell $(CXX) --version | egrep '^(Intel)'),)
ifneq ($(GPUCC),)
GPUFLAGS += -Wno-deprecated-builtins
endif
endif

# Avoid clang warning "overriding '-ffp-contract=fast' option with '-ffp-contract=on'" (#516)
# This patch does remove the warning, but I prefer to keep it disabled for the moment...
###ifneq ($(shell $(CXX) --version | egrep '^(clang|Apple clang|Intel)'),)
###$(BUILDDIR)/CrossSectionKernels_cpp.o: CXXFLAGS += -Wno-overriding-t-option
###ifneq ($(GPUCC),)
###$(BUILDDIR)/CrossSectionKernels_$(GPUSUFFIX).o: GPUFLAGS += $(XCOMPILERFLAG) -Wno-overriding-t-option
###endif
###endif

#### Apply special build flags only to CPPProcess.o (-flto)
###$(BUILDDIR)/CPPProcess_cpp.o: CXXFLAGS += -flto

#### Apply special build flags only to CPPProcess.o (AVXFLAGS)
###$(BUILDDIR)/CPPProcess_cpp.o: CXXFLAGS += $(AVXFLAGS)

# Generic target and build rules: objects from C++ compilation
# (NB do not include CUDA_INC here! add it only for NVTX or curand #679)
$(BUILDDIR)/%_cpp.o : %.cc *.h ../../src/*.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(INCFLAGS) $(CXXFLAGS) -c $< -o $@

# Generic target and build rules: objects from CUDA or HIP compilation
ifneq ($(GPUCC),)
$(BUILDDIR)/%_$(GPUSUFFIX).o : %.cc *.h ../../src/*.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(GPUCC) $(CPPFLAGS) $(INCFLAGS) $(GPUFLAGS) -c -x $(GPULANGUAGE) $< -o $@
endif

#-------------------------------------------------------------------------------

# Target (and build rules): common (src) library
commonlib : $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc $(BUILDDIR)/.build.$(TAG)
	$(MAKE) -C ../../src $(MAKEDEBUG) -f $(CUDACPP_SRC_MAKEFILE)

#-------------------------------------------------------------------------------

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')
###$(info processid_short=$(processid_short))

MG5AMC_CXXLIB = mg5amc_$(processid_short)_cpp
cxx_objects_lib=$(BUILDDIR)/CPPProcess_cpp.o $(BUILDDIR)/color_sum_cpp.o $(BUILDDIR)/MatrixElementKernels_cpp.o $(BUILDDIR)/BridgeKernels_cpp.o $(BUILDDIR)/CrossSectionKernels_cpp.o
cxx_objects_exe=$(BUILDDIR)/CommonRandomNumberKernel_cpp.o $(BUILDDIR)/RamboSamplingKernels_cpp.o

ifneq ($(GPUCC),)
MG5AMC_GPULIB = mg5amc_$(processid_short)_$(GPUSUFFIX)
gpu_objects_lib=$(BUILDDIR)/CPPProcess_$(GPUSUFFIX).o $(BUILDDIR)/color_sum_$(GPUSUFFIX).o $(BUILDDIR)/MatrixElementKernels_$(GPUSUFFIX).o $(BUILDDIR)/BridgeKernels_$(GPUSUFFIX).o $(BUILDDIR)/CrossSectionKernels_$(GPUSUFFIX).o
gpu_objects_exe=$(BUILDDIR)/CommonRandomNumberKernel_$(GPUSUFFIX).o $(BUILDDIR)/RamboSamplingKernels_$(GPUSUFFIX).o
endif

# Target (and build rules): C++ and CUDA/HIP shared libraries
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(BUILDDIR)/fbridge_cpp.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: cxx_objects_lib += $(BUILDDIR)/fbridge_cpp.o
$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib)
	$(CXX) -shared -o $@ $(cxx_objects_lib) $(CXXLIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

ifneq ($(GPUCC),)
$(LIBDIR)/lib$(MG5AMC_GPULIB).so: $(BUILDDIR)/fbridge_$(GPUSUFFIX).o
$(LIBDIR)/lib$(MG5AMC_GPULIB).so: gpu_objects_lib += $(BUILDDIR)/fbridge_$(GPUSUFFIX).o
$(LIBDIR)/lib$(MG5AMC_GPULIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(gpu_objects_lib)
	$(GPUCC) --shared -o $@ $(gpu_objects_lib) $(GPULIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB) $(BLASLIBFLAGS)
# Bypass std::filesystem completely to ease portability on LUMI #803
#ifneq ($(findstring hipcc,$(GPUCC)),)
#	$(GPUCC) --shared -o $@ $(gpu_objects_lib) $(GPULIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB) -lstdc++fs
#else
#	$(GPUCC) --shared -o $@ $(gpu_objects_lib) $(GPULIBFLAGSRPATH2) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)
#endif
endif

#-------------------------------------------------------------------------------

# Target (and build rules): Fortran include files
###$(INCDIR)/%.inc : ../%.inc
###	@if [ ! -d $(INCDIR) ]; then echo "mkdir -p $(INCDIR)"; mkdir -p $(INCDIR); fi
###	\cp $< $@

#-------------------------------------------------------------------------------

# Target (and build rules): C++ and CUDA/HIP standalone executables
###$(cxx_checkmain): LIBFLAGS += $(CXXLIBFLAGSASAN)
$(cxx_checkmain): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cxx_checkmain): $(BUILDDIR)/check_sa_cpp.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel_cpp.o $(BUILDDIR)/HiprandRandomNumberKernel_cpp.o
	$(CXX) -o $@ $(BUILDDIR)/check_sa_cpp.o $(OMPFLAGS) -ldl -pthread $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CXXLIB) $(cxx_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel_cpp.o $(BUILDDIR)/HiprandRandomNumberKernel_cpp.o $(RNDLIBFLAGS)

ifneq ($(GPUCC),)
###$(gpu_checkmain): LIBFLAGS += $(GPULIBFLAGSASAN)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(gpu_checkmain): LIBFLAGS += -lintlc # compile with icpx and link with GPUCC (undefined reference to `_intel_fast_memcpy')
$(gpu_checkmain): LIBFLAGS += -lsvml # compile with icpx and link with GPUCC (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(gpu_checkmain): LIBFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif
$(gpu_checkmain): LIBFLAGS += $(GPULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(gpu_checkmain): LIBFLAGS += $(BLASLIBFLAGS)
$(gpu_checkmain): $(BUILDDIR)/check_sa_$(GPUSUFFIX).o $(LIBDIR)/lib$(MG5AMC_GPULIB).so $(gpu_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel_$(GPUSUFFIX).o $(BUILDDIR)/HiprandRandomNumberKernel_$(GPUSUFFIX).o
	$(GPUCC) -o $@ $(BUILDDIR)/check_sa_$(GPUSUFFIX).o $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_GPULIB) $(gpu_objects_exe) $(BUILDDIR)/CurandRandomNumberKernel_$(GPUSUFFIX).o $(BUILDDIR)/HiprandRandomNumberKernel_$(GPUSUFFIX).o $(RNDLIBFLAGS)
endif

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from Fortran compilation
# (NB In this makefile, this only applies to fcheck_sa_fortran.o)
# (NB -fPIC was added to fix clang16 build #904, but this seems better for other cases too and is consistent to c++ and cuda builds)
$(BUILDDIR)/%_fortran.o : %.f *.inc
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(FC) -I. -fPIC -c $< -o $@

# Generic target and build rules: objects from Fortran compilation
###$(BUILDDIR)/%_fortran.o : %.f *.inc
###	@if [ ! -d $(INCDIR) ]; then echo "mkdir -p $(INCDIR)"; mkdir -p $(INCDIR); fi
###	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
###	$(FC) -I. -I$(INCDIR) -c $< -o $@

# Target (and build rules): Fortran standalone executables
###$(BUILDDIR)/fcheck_sa_fortran.o : $(INCDIR)/fbridge.inc

###$(cxx_fcheckmain): LIBFLAGS += $(CXXLIBFLAGSASAN)
ifeq ($(UNAME_S),Darwin)
$(cxx_fcheckmain): LIBFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(cxx_fcheckmain): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cxx_fcheckmain): $(BUILDDIR)/fcheck_sa_fortran.o $(BUILDDIR)/fsampler_cpp.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so $(cxx_objects_exe)
ifneq ($(findstring hipcc,$(GPUCC)),) # link fortran/c++/hip using $FC when hipcc is used #802
	$(FC) -o $@ $(BUILDDIR)/fcheck_sa_fortran.o $(OMPFLAGS) $(BUILDDIR)/fsampler_cpp.o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_CXXLIB) $(cxx_objects_exe) -lstdc++
else
	$(CXX) -o $@ $(BUILDDIR)/fcheck_sa_fortran.o $(OMPFLAGS) $(BUILDDIR)/fsampler_cpp.o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_CXXLIB) $(cxx_objects_exe)
endif

ifneq ($(GPUCC),)
###$(gpu_fcheckmain): LIBFLAGS += $(GPULIBFLAGSASAN)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(gpu_fcheckmain): LIBFLAGS += -lintlc # compile with icpx and link with GPUCC (undefined reference to `_intel_fast_memcpy')
$(gpu_fcheckmain): LIBFLAGS += -lsvml # compile with icpx and link with GPUCC (undefined reference to `__svml_cos4_l9')
endif
ifeq ($(UNAME_S),Darwin)
$(gpu_fcheckmain): LIBFLAGS += -L$(shell dirname $(shell $(FC) --print-file-name libgfortran.dylib)) # add path to libgfortran on Mac #375
endif
$(gpu_fcheckmain): LIBFLAGS += $(GPULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(gpu_fcheckmain): LIBFLAGS += $(BLASLIBFLAGS)
$(gpu_fcheckmain): $(BUILDDIR)/fcheck_sa_fortran.o $(BUILDDIR)/fsampler_$(GPUSUFFIX).o $(LIBDIR)/lib$(MG5AMC_GPULIB).so $(gpu_objects_exe)
ifneq ($(findstring hipcc,$(GPUCC)),) # link fortran/c++/hip using $FC when hipcc is used #802
	$(FC) -o $@ $(BUILDDIR)/fcheck_sa_fortran.o $(BUILDDIR)/fsampler_$(GPUSUFFIX).o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_GPULIB) $(gpu_objects_exe) -lstdc++ -L$(HIP_HOME)/lib -lamdhip64
else
	$(GPUCC) -o $@ $(BUILDDIR)/fcheck_sa_fortran.o $(BUILDDIR)/fsampler_$(GPUSUFFIX).o $(LIBFLAGS) -lgfortran -L$(LIBDIR) -l$(MG5AMC_GPULIB) $(gpu_objects_exe)
endif
endif

#-------------------------------------------------------------------------------

# Target (and build rules): test objects and test executable
ifeq ($(GPUCC),)
$(BUILDDIR)/testxxx_cpp.o: $(GTESTLIBS)
$(BUILDDIR)/testxxx_cpp.o: INCFLAGS += $(GTESTINC)
$(BUILDDIR)/testxxx_cpp.o: testxxx_cc_ref.txt
$(cxx_testmain): $(BUILDDIR)/testxxx_cpp.o
$(cxx_testmain): cxx_objects_exe += $(BUILDDIR)/testxxx_cpp.o # Comment out this line to skip the C++ test of xxx functions
else
$(BUILDDIR)/testxxx_$(GPUSUFFIX).o: $(GTESTLIBS)
$(BUILDDIR)/testxxx_$(GPUSUFFIX).o: INCFLAGS += $(GTESTINC)
$(BUILDDIR)/testxxx_$(GPUSUFFIX).o: testxxx_cc_ref.txt
$(gpu_testmain): $(BUILDDIR)/testxxx_$(GPUSUFFIX).o
$(gpu_testmain): gpu_objects_exe += $(BUILDDIR)/testxxx_$(GPUSUFFIX).o # Comment out this line to skip the CUDA/HIP test of xxx functions
endif

ifneq ($(UNAME_S),Darwin) # Disable testmisc on Darwin (workaround for issue #838)
ifeq ($(GPUCC),)
$(BUILDDIR)/testmisc_cpp.o: $(GTESTLIBS)
$(BUILDDIR)/testmisc_cpp.o: INCFLAGS += $(GTESTINC)
$(cxx_testmain): $(BUILDDIR)/testmisc_cpp.o
$(cxx_testmain): cxx_objects_exe += $(BUILDDIR)/testmisc_cpp.o # Comment out this line to skip the C++ miscellaneous tests
else
$(BUILDDIR)/testmisc_$(GPUSUFFIX).o: $(GTESTLIBS)
$(BUILDDIR)/testmisc_$(GPUSUFFIX).o: INCFLAGS += $(GTESTINC)
$(gpu_testmain): $(BUILDDIR)/testmisc_$(GPUSUFFIX).o
$(gpu_testmain): gpu_objects_exe += $(BUILDDIR)/testmisc_$(GPUSUFFIX).o # Comment out this line to skip the CUDA/HIP miscellaneous tests
endif
endif

ifeq ($(GPUCC),)
$(BUILDDIR)/runTest_cpp.o: $(GTESTLIBS)
$(BUILDDIR)/runTest_cpp.o: INCFLAGS += $(GTESTINC)
$(cxx_testmain): $(BUILDDIR)/runTest_cpp.o
$(cxx_testmain): cxx_objects_exe += $(BUILDDIR)/runTest_cpp.o
else
$(BUILDDIR)/runTest_$(GPUSUFFIX).o: $(GTESTLIBS)
$(BUILDDIR)/runTest_$(GPUSUFFIX).o: INCFLAGS += $(GTESTINC)
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(gpu_testmain): LIBFLAGS += -lintlc # compile with icpx and link with GPUCC (undefined reference to `_intel_fast_memcpy')
$(gpu_testmain): LIBFLAGS += -lsvml # compile with icpx and link with GPUCC (undefined reference to `__svml_cos4_l9')
else ifneq ($(shell $(CXX) --version | grep ^nvc++),) # support nvc++ #531
$(gpu_testmain): LIBFLAGS += -L$(patsubst %bin/nvc++,%lib,$(subst ccache ,,$(CXX))) -lnvhpcatm -lnvcpumath -lnvc
endif
$(gpu_testmain): $(BUILDDIR)/runTest_$(GPUSUFFIX).o
$(gpu_testmain): gpu_objects_exe  += $(BUILDDIR)/runTest_$(GPUSUFFIX).o
endif

ifeq ($(GPUCC),)
$(cxx_testmain): $(GTESTLIBS)
$(cxx_testmain): INCFLAGS +=  $(GTESTINC)
$(cxx_testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest # adding also -lgtest_main is no longer necessary since we added main() to testxxx.cc
else
$(gpu_testmain): $(GTESTLIBS)
$(gpu_testmain): INCFLAGS +=  $(GTESTINC)
$(gpu_testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest # adding also -lgtest_main is no longer necessary since we added main() to testxxx.cc
endif

ifeq ($(GPUCC),) # if at all, OMP is used only in CXX builds (not in GPU builds)
ifneq ($(OMPFLAGS),)
ifneq ($(shell $(CXX) --version | egrep '^Intel'),)
$(cxx_testmain): LIBFLAGS += -liomp5 # see #578 (not '-qopenmp -static-intel' as in https://stackoverflow.com/questions/45909648)
else ifneq ($(shell $(CXX) --version | egrep '^clang'),)
$(cxx_testmain): LIBFLAGS += -L $(shell dirname $(shell $(CXX) -print-file-name=libc++.so)) -lomp # see #604
###else ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
###$(cxx_testmain): LIBFLAGS += ???? # OMP is not supported yet by cudacpp for Apple clang (see #578 and #604)
else
$(cxx_testmain): LIBFLAGS += -lgomp
endif
endif
endif

# Test quadmath in testmisc.cc tests for constexpr_math #627
###ifeq ($(GPUCC),)
###$(cxx_testmain): LIBFLAGS += -lquadmath
###else
###$(gpu_testmain): LIBFLAGS += -lquadmath
###endif

# Bypass std::filesystem completely to ease portability on LUMI #803
###ifneq ($(findstring hipcc,$(GPUCC)),)
###$(gpu_testmain): LIBFLAGS += -lstdc++fs
###endif

ifeq ($(GPUCC),) # link only runTest_cpp.o
###$(cxx_testmain): LIBFLAGS += $(CXXLIBFLAGSASAN)
$(cxx_testmain): LIBFLAGS += $(CXXLIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(cxx_testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects_lib) $(cxx_objects_exe) $(GTESTLIBS)
	$(CXX) -o $@ $(cxx_objects_lib) $(cxx_objects_exe) -ldl -pthread $(LIBFLAGS)
else # link only runTest_$(GPUSUFFIX).o (new: in the past, this was linking both runTest_cpp.o and runTest_$(GPUSUFFIX).o)
###$(gpu_testmain): LIBFLAGS += $(GPULIBFLAGSASAN)
$(gpu_testmain): LIBFLAGS += $(GPULIBFLAGSRPATH) # avoid the need for LD_LIBRARY_PATH
$(gpu_testmain): LIBFLAGS += $(BLASLIBFLAGS)
$(gpu_testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(gpu_objects_lib) $(gpu_objects_exe) $(GTESTLIBS)
ifneq ($(findstring hipcc,$(GPUCC)),) # link fortran/c++/hip using $FC when hipcc is used #802
	$(FC) -o $@ $(gpu_objects_lib) $(gpu_objects_exe) -ldl $(LIBFLAGS) -lstdc++ -lpthread -L$(HIP_HOME)/lib -lamdhip64
else
	$(GPUCC) -o $@ $(gpu_objects_lib) $(gpu_objects_exe) -ldl $(LIBFLAGS) -lcuda
endif
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

# Target: build all targets in all BACKEND modes (each BACKEND mode in a separate build directory)
# Split the bldall target into separate targets to allow parallel 'make -j bldall' builds
# (Obsolete hack, no longer needed as there is no INCDIR: add a fbridge.inc dependency to bldall, to ensure it is only copied once for all BACKEND modes)
bldcuda:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cuda -f $(CUDACPP_MAKEFILE)

bldhip:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=hip -f $(CUDACPP_MAKEFILE)

bldnone:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cppnone -f $(CUDACPP_MAKEFILE)

bldsse4:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cppsse4 -f $(CUDACPP_MAKEFILE)

bldavx2:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cppavx2 -f $(CUDACPP_MAKEFILE)

bld512y:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cpp512y -f $(CUDACPP_MAKEFILE)

bld512z:
	@echo
	$(MAKE) USEBUILDDIR=1 BACKEND=cpp512z -f $(CUDACPP_MAKEFILE)

ifeq ($(UNAME_P),ppc64le)
###bldavxs: $(INCDIR)/fbridge.inc bldnone bldsse4
bldavxs: bldnone bldsse4
else ifeq ($(UNAME_P),arm)
###bldavxs: $(INCDIR)/fbridge.inc bldnone bldsse4
bldavxs: bldnone bldsse4
else
###bldavxs: $(INCDIR)/fbridge.inc bldnone bldsse4 bldavx2 bld512y bld512z
bldavxs: bldnone bldsse4 bldavx2 bld512y bld512z
endif

ifneq ($(HIP_HOME),)
ifneq ($(CUDA_HOME),)
bldall: bldhip bldcuda bldavxs
else
bldall: bldhip bldavxs
endif
else
ifneq ($(CUDA_HOME),)
bldall: bldcuda bldavxs
else
bldall: bldavxs
endif
endif

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
ifeq ($(USEBUILDDIR),1)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
	rm -f $(LIBDIR)/lib*.so
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
	@echo GPUCC=$(GPUCC)
ifneq ($(GPUCC),)
	$(GPUCC) --version
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

# Target: 'make test' (execute runTest.exe, and compare check.exe with fcheck.exe)
# [NB: THIS IS WHAT IS TESTED IN THE GITHUB CI!]
# [NB: This used to be called 'make check' but the name has been changed as this has nothing to do with 'check.exe']
test: runTest cmpFcheck

# Target: runTest (run the C++ or CUDA/HIP test executable runTest.exe)
runTest: all.$(TAG)
ifeq ($(GPUCC),)
	$(RUNTIME) $(BUILDDIR)/runTest_cpp.exe
else
	$(RUNTIME) $(BUILDDIR)/runTest_$(GPUSUFFIX).exe
endif

# Target: runCheck (run the C++ or CUDA/HIP standalone executable check.exe, with a small number of events)
runCheck: all.$(TAG)
ifeq ($(GPUCC),)
	$(RUNTIME) $(BUILDDIR)/check_cpp.exe -p 2 32 2
else
	$(RUNTIME) $(BUILDDIR)/check_$(GPUSUFFIX).exe -p 2 32 2
endif

# Target: runFcheck (run the Fortran standalone executable - with C++ or CUDA/HIP MEs - fcheck.exe, with a small number of events)
runFcheck: all.$(TAG)
ifeq ($(GPUCC),)
	$(RUNTIME) $(BUILDDIR)/fcheck_cpp.exe 2 32 2
else
	$(RUNTIME) $(BUILDDIR)/fcheck_$(GPUSUFFIX).exe 2 32 2
endif

# Target: cmpFcheck (compare ME results from the C++/CUDA/HIP and Fortran with C++/CUDA/HIP MEs standalone executables, with a small number of events)
cmpFcheck: all.$(TAG)
	@echo
ifeq ($(GPUCC),)
	@echo "$(BUILDDIR)/check_cpp.exe --common -p 2 32 2"
	@echo "$(BUILDDIR)/fcheck_cpp.exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(BUILDDIR)/check_cpp.exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(BUILDDIR)/fcheck_cpp.exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/C++)    = $${me1}"; echo "Avg ME (F77/C++)    = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/C++) returned NaN"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/C++) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi
else
	@echo "$(BUILDDIR)/check_$(GPUSUFFIX).exe --common -p 2 32 2"
	@echo "$(BUILDDIR)/fcheck_$(GPUSUFFIX).exe 2 32 2"
	@me1=$(shell $(RUNTIME) $(BUILDDIR)/check_$(GPUSUFFIX).exe --common -p 2 32 2 | grep MeanMatrix | awk '{print $$4}'); me2=$(shell $(RUNTIME) $(BUILDDIR)/fcheck_$(GPUSUFFIX).exe 2 32 2 | grep Average | awk '{print $$4}'); echo "Avg ME (C++/GPU)   = $${me1}"; echo "Avg ME (F77/GPU)   = $${me2}"; if [ "$${me2}" == "NaN" ]; then echo "ERROR! Fortran calculation (F77/GPU) crashed"; elif [ "$${me2}" == "" ]; then echo "ERROR! Fortran calculation (F77/GPU) crashed"; else python3 -c "me1=$${me1}; me2=$${me2}; reldif=abs((me2-me1)/me1); print('Relative difference =', reldif); ok = reldif <= 2E-4; print ( '%s (relative difference %s 2E-4)' % ( ('OK','<=') if ok else ('ERROR','>') ) ); import sys; sys.exit(0 if ok else 1)"; fi
endif

# Target: cuda-memcheck (run the CUDA standalone executable gcheck.exe with a small number of events through cuda-memcheck)
cuda-memcheck: all.$(TAG)
	$(RUNTIME) $(CUDA_HOME)/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes $(BUILDDIR)/check_$(GPUSUFFIX).exe -p 2 32 2

#-------------------------------------------------------------------------------
