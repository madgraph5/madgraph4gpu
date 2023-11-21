# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: J. Teig, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

#=== Determine the name of this makefile (https://ftp.gnu.org/old-gnu/Manuals/make-3.80/html_node/make_17.html)
#=== NB: assume that the same name (e.g. cudacpp.mk, Makefile...) is used in the Subprocess and src directories

THISMK = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))

#-------------------------------------------------------------------------------

#=== Use bash in the Makefile (https://www.gnu.org/software/make/manual/html_node/Choosing-the-Shell.html)

SHELL := /bin/bash

#-------------------------------------------------------------------------------

#=== Configure the C++ compiler

include ../Source/make_opts

MG_CXXFLAGS += -fPIC -I. $(USE_NVTX)
ifeq ($(shell $(CXX) --version | grep ^nvc++),)
MG_CXXFLAGS += -ffast-math # see issue #117
endif

# Note: AR, CXX and FC are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
###RANLIB = ranlib

#-------------------------------------------------------------------------------

#=== Configure the CUDA compiler (note: NVCC is already exported including ccache)

###$(info NVCC=$(NVCC))

#-------------------------------------------------------------------------------

#=== Configure ccache for C++ builds (note: NVCC is already exported including ccache)

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif
#ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
#  override AR:=ccache $(AR)
#endif

#-------------------------------------------------------------------------------

#=== Configure PowerPC-specific compiler flags for CUDA and C++

# Assuming uname is available, detect if architecture is PowerPC
UNAME_P := $(shell uname -p)

# PowerPC-specific CXX compiler flags (being reviewed)
ifeq ($(UNAME_P),ppc64le)
  MG_CXXFLAGS+= -mcpu=power9 -mtune=power9 # gains ~2-3% both for none and sse4
  # Throughput references without the extra flags below: none=1.41-1.42E6, sse4=2.15-2.19E6
  ###MG_CXXFLAGS+= -DNO_WARN_X86_INTRINSICS # no change
  ###MG_CXXFLAGS+= -fpeel-loops # no change
  ###MG_CXXFLAGS+= -funroll-loops # gains ~1% for none, loses ~1% for sse4
  ###MG_CXXFLAGS+= -ftree-vectorize # no change
  ###MG_CXXFLAGS+= -flto # BUILD ERROR IF THIS ADDED IN SRC?!
else
  ###AR=gcc-ar # needed by -flto
  ###RANLIB=gcc-ranlib # needed by -flto
  ###MG_CXXFLAGS+= -flto # NB: build error from src/Makefile unless gcc-ar and gcc-ranlib are used
  ######MG_CXXFLAGS+= -fno-semantic-interposition # no benefit (neither alone, nor combined with -flto)
endif

#-------------------------------------------------------------------------------

#=== Set the CUDA/C++ compiler flags appropriate to user-defined choices of AVX, FPTYPE, HELINL, HRDCOD, RNDGEN

# Set the build flags appropriate to OMPFLAGS
###$(info OMPFLAGS=$(OMPFLAGS))
MG_CXXFLAGS += $(OMPFLAGS)

# Set the build flags appropriate to each AVX choice (example: "make AVX=none")
# [NB MGONGPU_PVW512 is needed because "-mprefer-vector-width=256" is not exposed in a macro]
# [See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96476]
ifeq ($(NVCC),)
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
    else ifneq ($(AVX),none)
      $(error Unknown AVX='$(AVX)': only 'none', 'sse4', 'avx2', '512y' and '512z' are supported)
    endif
  endif
  # For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations?
  MG_CXXFLAGS+= $(AVXFLAGS)
endif

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
###$(info FPTYPE=$(FPTYPE))
ifeq ($(FPTYPE),d)
  MG_CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_DOUBLE
else ifeq ($(FPTYPE),f)
  MG_CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT -DMGONGPU_FPTYPE2_FLOAT
else ifeq ($(FPTYPE),m)
  MG_CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE -DMGONGPU_FPTYPE2_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'd', 'f' and 'm' are supported)
endif

# Set the build flags appropriate to each HELINL choice (example: "make HELINL=1")
###$(info HELINL=$(HELINL))
ifeq ($(HELINL),1)
  MG_CXXFLAGS += -DMGONGPU_INLINE_HELAMPS
else ifneq ($(HELINL),0)
  $(error Unknown HELINL='$(HELINL)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each HRDCOD choice (example: "make HRDCOD=1")
###$(info HRDCOD=$(HRDCOD))
ifeq ($(HRDCOD),1)
  MG_CXXFLAGS += -DMGONGPU_HARDCODE_PARAM
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each RNDGEN choice (example: "make RNDGEN=hasNoCurand")
###$(info RNDGEN=$(RNDGEN))
ifeq ($(RNDGEN),hasNoCurand)
  MG_CXXFLAGS += -DMGONGPU_HAS_NO_CURAND
else ifneq ($(RNDGEN),hasCurand)
  $(error Unknown RNDGEN='$(RNDGEN)': only 'hasCurand' and 'hasNoCurand' are supported)
endif

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
ifneq ($(NVCC),)
  override DIRTAG = cuda_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)
else
  override DIRTAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)
endif

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of CUDA and no-CUDA environment builds with different random number generators)
ifneq ($(NVCC),)
  override TAG = cuda_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)_$(RNDGEN)
else
  override TAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)_$(RNDGEN)
endif

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
###$(info Current directory is $(shell pwd))
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
  override LIBDIRREL = ../lib/$(BUILDDIR)
  ###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR=1 is set))
else
  override BUILDDIR = .
  override LIBDIRREL = ../lib
  ###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is not set))
endif
######$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

# Workaround for Mac #375 (I did not manage to fix rpath with @executable_path): use absolute paths for LIBDIR
# (NB: this is quite ugly because it creates the directory if it does not exist - to avoid removing src by mistake)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
override LIBDIR = $(shell mkdir -p $(LIBDIRREL); cd $(LIBDIRREL); pwd)
ifeq ($(wildcard $(LIBDIR)),)
$(error Directory LIBDIR="$(LIBDIR)" should have been created by now)
endif
else
override LIBDIR = $(LIBDIRREL)
endif

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

# NB1: there are no CUDA targets in src as we avoid RDC!
# NB2: CUDA includes for curand.h are no longer needed in the C++ code anywhere in src!

MG5AMC_COMMONLIB = mg5amc_common

# First target (default goal)
all.$(TAG): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from C++ compilation
$(BUILDDIR)/%.o : %.cc *.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(MG_CXXFLAGS) $(CXXFLAGS) -c $< -o $@

# Generic target and build rules: objects from CUDA compilation
$(BUILDDIR)/%_cu.o : %.cc *.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(MG_NVCCFLAGS) $(NVCCFLAGS) -c -x cu $< -o $@

#-------------------------------------------------------------------------------

cxx_objects=$(addprefix $(BUILDDIR)/, Parameters_sm.o read_slha.o)
ifneq ($(NVCC),)
cu_objects=$(addprefix $(BUILDDIR)/, Parameters_sm_cu.o)
endif

# Target (and build rules): common (src) library
ifneq ($(NVCC),)
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects) $(cu_objects)
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	$(NVCC) -shared -o $@ $(cxx_objects) $(cu_objects)
else
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects)
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	$(CXX) $(MG_LDFLAGS) $(LDFLAGS) -shared -o $@ $(cxx_objects)
endif

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

BUILD_DIRS := $(wildcard build.*)
NUM_BUILD_DIRS := $(words $(BUILD_DIRS))

clean:
ifeq ($(USEBUILDDIR),1)
ifeq ($(NUM_BUILD_DIRS),1)
	$(info USEBUILDDIR=1, only one src build directory found.)
	rm -rf ../lib/$(BUILD_DIRS)
	rm -rf $(BUILD_DIRS)
else ifeq ($(NUM_BUILD_DIRS),0)
	$(error USEBUILDDIR=1, but no src build directories are found.)
else
	$(error Multiple src BUILDDIR's found! Use 'cleannone', 'cleansse4', 'cleanavx2', 'clean512y','clean512z', 'cleancuda' or 'cleanall'.)
endif
else
	rm -f ../lib/lib$(MG5AMC_COMMONLIB).so
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

cleanall:
	@echo
	rm -f ../lib/lib$(MG5AMC_COMMONLIB).so
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
	@echo
	rm -rf ../lib/build.*
	rm -rf build.*

# Target: clean different builds

cleannone:
	rm -rf ../lib/build.none_*
	rm -rf build.none_*

cleansse4:
	rm -rf ../lib/build.sse4_*
	rm -rf build.sse4_*

cleanavx2:
	rm -rf ../lib/build.avx2_*
	rm -rf build.avx2_*

clean512y:
	rm -rf ../lib/build.512y_*
	rm -rf build.512y_*

clean512z:
	rm -rf ../lib/build.512z_*
	rm -rf build.512z_*

cleancuda:
	rm -rf ../lib/build.cuda_*
	rm -rf build.cuda_*

cleandir:
	rm -f ./*.o ./*.exe
	rm -f ../lib/lib$(MG5AMC_COMMONLIB).so

#-------------------------------------------------------------------------------
