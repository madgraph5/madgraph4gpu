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

# Add option to enable CI profiler use
$(info ENABLE_CI_PROFILER=$(ENABLE_CI_PROFILER))
ifeq ($(ENABLE_CI_PROFILER),1)
  CXXFLAGS += --gcc-toolchain="/cvmfs/sft.cern.ch/lcg/releases/gcc/11.3.0-ad0f5/x86_64-centos8"
endif

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
$(BUILDDIR)/%.o : %.cc *.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------

cxx_objects=$(addprefix $(BUILDDIR)/, Parameters_sm.o read_slha.o)

# Target (and build rules): common (src) library
$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so : $(cxx_objects)
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(CXX) -shared -o$@ $(cxx_objects)

#-------------------------------------------------------------------------------

# Target (and build rules): install libraries and headers (for use by MadEvent in Fortran)
INSTALL_HEADERS=Parameters_sm.h mgOnGpuConfig.h mgOnGpuFptypes.h mgOnGpuCxtypes.h mgOnGpuVectors.h read_slha.h extras.h
INSTALL_INC_DIR=../include

install: all.$(TAG) $(INSTALL_INC_DIR) $(addprefix $(INSTALL_INC_DIR)/, $(INSTALL_HEADERS))

$(INSTALL_INC_DIR) :
	if [ ! -d $(INSTALL_INC_DIR) ]; then mkdir $(INSTALL_INC_DIR); fi

$(INSTALL_INC_DIR)/%.h : %.h
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
