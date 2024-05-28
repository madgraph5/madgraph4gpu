# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
# Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi (2020-2024) for the MG5aMC CUDACPP plugin.

#=== Determine the name of this makefile (https://ftp.gnu.org/old-gnu/Manuals/make-3.80/html_node/make_17.html)
#=== NB: assume that the same name (e.g. cudacpp.mk, Makefile...) is used in the Subprocess and src directories

THISMK = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))

#-------------------------------------------------------------------------------

#=== Use bash in the Makefile (https://www.gnu.org/software/make/manual/html_node/Choosing-the-Shell.html)

SHELL := /bin/bash

#-------------------------------------------------------------------------------

#=== Configure common compiler flags for CUDA and C++

INCFLAGS = -I.

#-------------------------------------------------------------------------------

#=== Configure the C++ compiler (note: CXXFLAGS has been exported from cudacpp.mk)

###$(info CXXFLAGS=$(CXXFLAGS))

# Note: AR, CXX and FC are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
###RANLIB = ranlib

# Add -mmacosx-version-min=11.3 to avoid "ld: warning: object file was built for newer macOS version than being linked"
LDFLAGS =
ifneq ($(shell $(CXX) --version | egrep '^Apple clang'),)
  LDFLAGS += -mmacosx-version-min=11.3
endif

#-------------------------------------------------------------------------------

#=== Configure the GPU (CUDA or HIP) compiler (note: GPUCC including ccache, GPUFLAGS, GPULANGUAGE, GPUSUFFIX have been exported from cudacpp.mk)

###$(info GPUCC=$(GPUCC))
###$(info GPUFLAGS=$(GPUFLAGS))
###$(info GPULANGUAGE=$(GPULANGUAGE))
###$(info GPUSUFFIX=$(GPUSUFFIX))

#-------------------------------------------------------------------------------

#=== Configure ccache for C++ builds (note: GPUCC has been exported from cudacpp.mk including ccache)

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif
#ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
#  override AR:=ccache $(AR)
#endif

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Use the build directory exported from cudacpp.mk
###$(info CUDACPP_BUILDDIR=$(CUDACPP_BUILDDIR))

# Use the build lockfile "full" tag exported from cudacpp.mk
###$(info TAG=$(TAG))

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
###$(info Current directory is $(shell pwd))
override BUILDDIR = $(CUDACPP_BUILDDIR)
ifeq ($(USEBUILDDIR),1)
  override LIBDIRREL = ../lib/$(BUILDDIR)
  ###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR=1 is set))
else
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

ifeq ($(GPUCC),)
MG5AMC_COMMONLIB = mg5amc_common_cpp
else
MG5AMC_COMMONLIB = mg5amc_common_$(GPUSUFFIX)
endif

# Explicitly define the default goal (this is not necessary as it is the first target, which is implicitly the default goal)
.DEFAULT_GOAL := all.$(TAG)

# First target (default goal)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/.build.$(TAG) $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

# Target (and build options): debug
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
override oldtagsl=`if [ -d $(LIBDIR) ]; then find $(LIBDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`

$(BUILDDIR)/.build.$(TAG): $(LIBDIR)/.build.$(TAG)

$(LIBDIR)/.build.$(TAG):
	@if [ "$(oldtagsl)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist in $(LIBDIR) for other tags:\n$(oldtagsl)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make cleanall'"; exit 1; fi
	@if [ "$(oldtagsb)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist in $(BUILDDIR) for other tags:\n$(oldtagsb)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make cleanall'"; exit 1; fi
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	@touch $(LIBDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@touch $(BUILDDIR)/.build.$(TAG)

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from C++ compilation
$(BUILDDIR)/%_cpp.o : %.cc *.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(INCFLAGS) $(CXXFLAGS) -c $< -o $@

# Generic target and build rules: objects from CUDA compilation
ifneq ($(GPUCC),)
$(BUILDDIR)/%_$(GPUSUFFIX).o : %.cc *.h $(BUILDDIR)/.build.$(TAG)
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(GPUCC) $(CPPFLAGS) $(INCFLAGS) $(GPUFLAGS) -c -x $(GPULANGUAGE) $< -o $@
endif

#-------------------------------------------------------------------------------

cxx_objects=$(addprefix $(BUILDDIR)/, read_slha_cpp.o)
ifeq ($(GPUCC),)
  cxx_objects+=$(addprefix $(BUILDDIR)/, Parameters_sm_cpp.o)
else
  gpu_objects=$(addprefix $(BUILDDIR)/, Parameters_sm_$(GPUSUFFIX).o)
endif

# Target (and build rules): common (src) library
ifeq ($(GPUCC),)
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects)
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	$(CXX) -shared -o $@ $(cxx_objects) $(LDFLAGS)
else
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects) $(gpu_objects)
	@if [ ! -d $(LIBDIR) ]; then echo "mkdir -p $(LIBDIR)"; mkdir -p $(LIBDIR); fi
	$(GPUCC) -shared -o $@ $(cxx_objects) $(gpu_objects) $(LDFLAGS)
endif

#-------------------------------------------------------------------------------

# Target: clean the builds
.PHONY: clean

clean:
ifeq ($(USEBUILDDIR),1)
	rm -rf $(LIBDIR)
	rm -rf $(BUILDDIR)
else
	rm -f $(LIBDIR)/.build.* $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

cleanall:
	@echo
	$(MAKE) clean -f $(THISMK)
	@echo
	rm -rf $(LIBDIR)/build.*
	rm -rf build.*

#-------------------------------------------------------------------------------
