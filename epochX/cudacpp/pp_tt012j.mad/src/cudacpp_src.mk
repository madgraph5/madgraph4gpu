# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: J. Teig, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

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

#=== Set the CUDA/C++ compiler flags appropriate to user-defined choices of AVX, FPTYPE, HELINL, HRDCOD, RNDGEN

# Set the build flags appropriate to OMPFLAGS
###$(info OMPFLAGS=$(OMPFLAGS))
MG_CXXFLAGS += $(OMPFLAGS)

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
DIRTAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of CUDA and no-CUDA environment builds with different random number generators)
TAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)_$(RNDGEN)

# Build directory:
BUILDDIR := build.$(DIRTAG)
LIBDIRREL := ../lib/$(BUILDDIR)

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
ifeq ($(AVX),cuda)
COMPILER=$(NVCC)
cu_objects=$(addprefix $(BUILDDIR)/, Parameters_sm_cu.o)
else
COMPILER=$(CXX)
cu_objects=
endif

# Target (and build rules): common (src) library
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects) $(cu_objects)
	mkdir -p $(LIBDIR)
	$(COMPILER) -shared -o $@ $(cxx_objects) $(cu_objects) $(MG_LDFLAGS) $(LDFLAGS)

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
	$(RM) -f ../lib/build.*/*.so
	$(RM) -rf build.*

#-------------------------------------------------------------------------------
