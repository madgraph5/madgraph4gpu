# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
#
# Copyright (C) 2021-2023 Argonne National Laboratory.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Modified by: N. Nichols for the MG5aMC SYCL plugin.

#=== Configure common compiler flags for SYCL build

INCFLAGS = -I.
OPTFLAGS = -O3 -march=native # this ends up in 

#-------------------------------------------------------------------------------

#=== Configure the SYCL compiler

CXXFLAGS = $(OPTFLAGS) -std=c++20 $(INCFLAGS) -fPIC -Wall -Wshadow -Wextra
CXXFLAGS+= -ffast-math # see issue #117

# Note: AR and CXX are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

#-------------------------------------------------------------------------------

#=== Configure ccache for SYCL build

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif

#-------------------------------------------------------------------------------

#=== Set the SYCL compiler flags appropriate to user-defined choices of FPTYPE, HELINL, HRDCOD

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'f' and 'd' are supported)
endif

# Set the build flags appropriate to each HELINL choice (example: "make HELINL=1")
ifeq ($(HELINL),1)
  CXXFLAGS += -DMGONGPU_INLINE_HELAMPS
else ifneq ($(HELINL),0)
  $(error Unknown HELINL='$(HELINL)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each HRDCOD choice (example: "make HRDCOD=1")
ifeq ($(HRDCOD),1)
  CXXFLAGS += -DMGONGPU_HARDCODE_PARAM
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each NTPBMAX choice (example: "make NTPBMAX=1024")
$(info NTPBMAX=$(NTPBMAX))
CXXFLAGS += -DMGONGPU_NTPBMAX=$(NTPBMAX)

# Set the build flags appropriate to each VECLVL choice (example: "make VECLVL=1")
$(info VECLVL=$(VECLVL))
ifeq ($(VECLVL),1)
  CXXFLAGS += -DMGONGPU_VEC_DIM=1
else ifeq ($(VECLVL),2)
  CXXFLAGS += -DMGONGPU_VEC_DIM=2
else ifeq ($(VECLVL),4)
  CXXFLAGS += -DMGONGPU_VEC_DIM=4
else ifeq ($(VECLVL),8)
  CXXFLAGS += -DMGONGPU_VEC_DIM=8
else ifeq ($(VECLVL),16)
  CXXFLAGS += -DMGONGPU_VEC_DIM=16
else
  $(error Unknown VECLVL='$(VECLVL)': only '1', '2', '4', '8', and '16' are supported)
endif

# Set the build flags appropriate to each CXTYPE choice (example: "make CXTYPE=std")
$(info CXTYPE=$(CXTYPE))
ifeq ($(CXTYPE),smpl)
  CXXFLAGS += -DMGONGPU_COMPLEX_CXSMPL
else ifeq ($(CXTYPE),extras)
  CXXFLAGS += -DMGONGPU_COMPLEX_EXTRAS
else ifeq ($(CXTYPE),std)
  CXXFLAGS += -DMGONGPU_COMPLEX_STD
else ifeq ($(CXTYPE),oneapi)
  CXXFLAGS += -DMGONGPU_COMPLEX_ONEAPI
else ifeq ($(CXTYPE),thrust)
  CXXFLAGS += -DMGONGPU_COMPLEX_CUTHRUST
else ifeq ($(CXTYPE),cucomplex)
  CXXFLAGS += -DMGONGPU_COMPLEX_CUCOMPLEX
else ifeq ($(CXTYPE),syclcplx)
  CXXFLAGS += -DMGONGPU_COMPLEX_SYCLCPLX
else
  $(error Unknown CXTYPE='$(CXTYPE)': only 'smpl', 'extras', 'std', 'oneapi', 'thrust', 'cucomplex', and 'syclcplx' are supported)
endif

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
override DIRTAG = $(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
override TAG = $(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
  override LIBDIR   = ../lib/$(BUILDDIR)
  ###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR=1 is set))
else
  override BUILDDIR = .
  override LIBDIR   = ../lib
  ###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is not set))
endif
######$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

MG5AMC_COMMONLIB = mg5amc_common

# First target (default goal)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

# Target (and build options): debug
debug: OPTFLAGS = -g -O0 -DDEBUG2
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
override oldtagsl=`if [ -d $(LIBDIR) ]; then find $(LIBDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`

$(BUILDDIR)/.build.$(TAG): $(LIBDIR)/.build.$(TAG)

$(LIBDIR)/.build.$(TAG):
	@if [ "$(oldtagsl)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist in $(LIBDIR) for other tags:\n$(oldtagsl)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@if [ "$(oldtagsb)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist in $(BUILDDIR) for other tags:\n$(oldtagsb)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	@touch $(LIBDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@touch $(BUILDDIR)/.build.$(TAG)

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from C++ compilation
$(BUILDDIR)/%%.o : %%.cc *.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------

cxx_objects=$(addprefix $(BUILDDIR)/, Parameters_%(model)s.o read_slha.o)

# Target (and build rules): common (src) library
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects)
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(CXX) -shared -o$@ $(cxx_objects)

#-------------------------------------------------------------------------------

# Target (and build rules): install libraries and headers (for use by MadEvent in Fortran)
INSTALL_HEADERS=Parameters_%(model)s.h mgOnGpuConfig.h mgOnGpuFptypes.h mgOnGpuCxtypes.h mgOnGpuVectors.h read_slha.h extras.h sycl_ext_complex.hpp
INSTALL_INC_DIR=../include

install: all.$(TAG) $(INSTALL_INC_DIR) $(addprefix $(INSTALL_INC_DIR)/, $(INSTALL_HEADERS))

$(INSTALL_INC_DIR) :
	if [ ! -d $(INSTALL_INC_DIR) ]; then mkdir $(INSTALL_INC_DIR); fi

$(INSTALL_INC_DIR)/%%.h : %%.h
	@if [ ! -d $(INSTALL_INC_DIR) ]; then mkdir $(INSTALL_INC_DIR); fi
	cp $< $@

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
	rm -rf $(LIBDIR)
ifneq ($(BUILDDIR),.)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif
	if [ -d $(INSTALL_INC_DIR) ]; then rm -rf $(INSTALL_INC_DIR); fi

cleanall:
	@echo
	make clean
	@echo
	rm -rf build.*

#-------------------------------------------------------------------------------
