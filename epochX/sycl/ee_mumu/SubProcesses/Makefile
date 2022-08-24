#=== Configure common compiler flags for CUDA and C++

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

# Note: AR and CXX are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

#-------------------------------------------------------------------------------

#=== Configure ccache for SYCL builds

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
  override CXX:=ccache $(CXX)
endif

#-------------------------------------------------------------------------------

#=== Configure defaults and check if user-defined choices exist for PTYPE, HELINL

# Set the default FPTYPE (floating point type) choice
ifeq ($(FPTYPE),)
  override FPTYPE = d
endif

# Set the default HELINL (inline helicities?) choice
ifeq ($(HELINL),)
  override HELINL = 0
endif

# Set the default NTPBMAX (maximum threads per block) choice
ifeq ($(NTPBMAX),)
  override NTPBMAX = 1024
endif

# Export FPTYPE, HELINL so that it is not necessary to pass them to the src Makefile too
export FPTYPE
export HELINL
export NTPBMAX

#-------------------------------------------------------------------------------

#=== Set the SYCL compiler flags appropriate to user-defined choices of FPTYPE, HELINL

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

# Set the build flags appropriate to each NTPBMAX choice (example: "make NTPBMAX=1024")
$(info NTPBMAX=$(NTPBMAX))
CXXFLAGS += -DMGONGPU_NTPBMAX=$(NTPBMAX)

#-------------------------------------------------------------------------------

#=== Configure build directories and build lockfiles ===

# Build directory "short" tag (defines target and path to the optional build directory)
override DIRTAG = $(FPTYPE)_inl$(HELINL)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
override TAG = $(FPTYPE)_inl$(HELINL)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
  override LIBDIR   = ../../lib/$(BUILDDIR)
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is set = 1))
else
  override BUILDDIR = .
  override LIBDIR   = ../../lib
  $(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG) (USEBUILDDIR is not set))
endif
###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

#===============================================================================
#=== Makefile TARGETS and build rules below
#===============================================================================

cu_main=$(BUILDDIR)/check.exe

# First target (default goal)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) ../../src/$(BUILDDIR)/.build.$(TAG) $(cu_main)

# Target (and build options): debug
MAKEDEBUG=
debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: MAKEDEBUG := debug
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtagsb)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtagsb)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

../../src/$(BUILDDIR)/.build.$(TAG):
	$(MAKE) -C ../../src $(MAKEDEBUG)

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from SYCL compilation
$(BUILDDIR)/%.o : %.cu *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -c $< -o $@

$(BUILDDIR)/%_cu.o : %.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -c $< -o $@

# Generic target and build rules: objects from SYCL compilation
$(BUILDDIR)/%.o : %.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -fPIC -c $< -o $@

# Apply special build flags only to CrossSectionKernel.cc (no fast math, see #117)
$(BUILDDIR)/CrossSectionKernels.o: CXXFLAGS += -fno-fast-math

#-------------------------------------------------------------------------------

# Target (and build rules): common (src) library
commonlib : $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc
	$(MAKE) -C ../../src $(MAKEDEBUG)

#-------------------------------------------------------------------------------

MG5AMC_CULIB = mg5amc_cu
cu_objects=$(BUILDDIR)/CPPProcess.o $(BUILDDIR)/RandomNumberKernels.o $(BUILDDIR)/RamboSamplingKernels.o $(BUILDDIR)/MatrixElementKernels.o $(BUILDDIR)/CrossSectionKernels.o

# Target (and build rules): SYCL shared libraries
$(LIBDIR)/lib$(MG5AMC_CULIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_objects)
	$(CXX) $(SYCLFLAGS) -shared -o $@ $(cu_objects) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

#-------------------------------------------------------------------------------

# Target (and build rules): SYCL standalone executables
$(cu_main): $(LIBDIR)/lib$(MG5AMC_CULIB).so
	$(CXX) $(CXXFLAGS) $(SYCLFLAGS) -o $@ check_sa.cc -pthread $(LIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CULIB) -lstdc++fs

#-------------------------------------------------------------------------------

# Target (and build rules): install libraries and headers (for use by MadEvent in Fortran)
INSTALL_HEADERS=CPPProcess.h MatrixElementKernels.h MemoryBuffers.h extras.h
INSTALL_INC_DIR=../../include
INSTALL_OBJECTS=CPPProcess.o MatrixElementKernels.o

install: all.$(TAG) $(INSTALL_INC_DIR) $(addprefix $(INSTALL_INC_DIR)/, $(INSTALL_HEADERS)) $(addprefix $(LIBDIR)/, $(INSTALL_OBJECTS))

$(INSTALL_INC_DIR) :
	$(MAKE) -C ../../src install

$(INSTALL_INC_DIR)/%.h : %.h
	@if [ ! -d $(INSTALL_INC_DIR) ]; then mkdir $(INSTALL_INC_DIR); fi
	cp $< $@

$(LIBDIR)/%.o : $(BUILDDIR)/%.o
	cp $< $@

#-------------------------------------------------------------------------------


# Target: clean the builds
.PHONY: clean

clean:
	make -C ../../src clean
ifneq ($(BUILDDIR),.)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

cleanall:
	@echo
	make USEBUILDDIR=0 clean
	@echo
	make USEBUILDDIR=0 -C ../../src cleanall
	rm -rf build.*

#-------------------------------------------------------------------------------

# Target: show system and compiler information
info:
	@echo ""
	@echo -n "hostname="
	@hostname
	@cat /proc/cpuinfo | grep "model name" | sort -u
	@cat /proc/cpuinfo | grep "flags" | sort -u
	@cat /proc/cpuinfo | grep "cpu cores" | sort -u
	@cat /proc/cpuinfo | grep "physical id" | sort -u
	@echo ""
	@echo USECCACHE=$(USECCACHE)
ifeq ($(USECCACHE),1)
	ccache --version | head -1
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

#-------------------------------------------------------------------------------

# Target: check (run the C++ test executable)
# [NB THIS IS WHAT IS USED IN THE GITHUB CI!]
check: runGcheck

# Target: runGcheck (run the SYCL standalone executable check.exe with a small number of events)
runGcheck: $(BUILDDIR)/check.exe
	DYLD_LIBRARY_PATH="$(LIBDIR):$(DYLD_LIBRARY_PATH)" LD_LIBRARY_PATH="$(LIBDIR):$(LD_LIBRARY_PATH)" $(BUILDDIR)/check.exe -p --device_id 0 2 32 2

#-------------------------------------------------------------------------------
