#=== Configure common compiler flags for CUDA and C++

INCFLAGS = -I.
OPTFLAGS = -O3 # this ends up in CUFLAGS too (should it?), cannot add -Ofast or -ffast-math here

# Dependency on src directory
MG5AMC_COMMONLIB = mg5amc_common
LIBFLAGS = -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)
INCFLAGS += -I../../src

# Dependency on tools directory
TOOLSDIR = ../../../../../tools
INCFLAGS += -I$(TOOLSDIR)

# Dependency on test directory
TESTDIR  = ../../../../../test
GTESTLIBDIR = $(TESTDIR)/googletest/build/lib/
GTESTLIBS   = $(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a

#-------------------------------------------------------------------------------

#=== Configure the C++ compiler

CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
###CXXFLAGS+= -g # FOR DEBUGGING ONLY

# Note: AR and CXX are implicitly defined if not set externally
# See https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

#-------------------------------------------------------------------------------

#=== Configure the CUDA compiler

# If CUDA_HOME is not set, try to set it from the location of nvcc
ifndef CUDA_HOME
  CUDA_HOME = $(patsubst %%bin/nvcc,%%,$(shell which nvcc 2>/dev/null))
  $(warning CUDA_HOME was not set: using "$(CUDA_HOME)")
endif

# Set NVCC as $(CUDA_HOME)/bin/nvcc if it exists
ifneq ($(wildcard $(CUDA_HOME)/bin/nvcc),)
  NVCC = $(CUDA_HOME)/bin/nvcc
  CUARCHNUM=70 # For V100 (Volta), e.g. SCALE cluster and CERN lxbatch, see https://developer.nvidia.com/cuda-gpus#compute
  ###CUARCHNUM=35 # For Jetson Nano (Maxwell), see https://developer.nvidia.com/cuda-gpus#compute
  USE_NVTX ?=-DUSE_NVTX
  CUARCHFLAGS = -arch=compute_$(CUARCHNUM)
  ###CUARCHFLAGS = -gencode arch=compute_$(CUARCHNUM),code=sm_$(CUARCHNUM)
  CUINC       = -I$(CUDA_HOME)/include/
  CULIBFLAGS  = -L$(CUDA_HOME)/lib64/ -lcurand # NB: -lcuda is not needed here!
  CUOPTFLAGS  = -lineinfo
  CUFLAGS     = $(OPTFLAGS) $(CUOPTFLAGS) $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math
  ###CUFLAGS    += -Xcompiler -Wall -Xcompiler -Wextra -Xcompiler -Wshadow
  ###NVCC_VERSION = $(shell $(NVCC) --version | grep 'Cuda compilation tools' | cut -d' ' -f5 | cut -d, -f1)
  CUFLAGS += -std=c++17 # need CUDA >= 11.2 (see #333): this is enforced in mgOnGpuConfig.h
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)
else
  # No cuda. Switch cuda compilation off and go to common random numbers in C++
  $(warning CUDA_HOME is not set or is invalid. Export CUDA_HOME to compile with cuda)
  override NVCC=
  override USE_NVTX=
  override CULIBFLAGS=
endif

# Set the host C++ compiler for nvcc
CUFLAGS += -ccbin $(shell which $(subst ccache ,,$(CXX)))

#-------------------------------------------------------------------------------

#=== Configure ccache for CUDA and C++ builds

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

#=== Configure PowerPC-specific compiler flags for CUDA and C++

# Assuming uname is available, detect if architecture is PowerPC
UNAME_P := $(shell uname -p)

# PowerPC-specific CXX compiler flags (being reviewed)
ifeq ($(UNAME_P),ppc64le)
  CXXFLAGS+= -mcpu=power9 -mtune=power9 # gains ~2-3%% both for none and sse4
  # Throughput references without the extra flags below: none=1.41-1.42E6, sse4=2.15-2.19E6
  ###CXXFLAGS+= -DNO_WARN_X86_INTRINSICS # no change
  ###CXXFLAGS+= -fpeel-loops # no change
  ###CXXFLAGS+= -funroll-loops # gains ~1%% for none, loses ~1%% for sse4
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
ifneq ($(shell $(CXX) --version | grep ^Intel),)
override OMPFLAGS = # disable OpenMP on the Intel compiler (on gcc this requires gcc>=9.3, issue #269)
endif
###OMPFLAGS ?= -fopenmp # TEMPORARELY DISABLE OMP (need to reassess MT)
override OMPFLAGS = # TEMPORARELY DISABLE OMP (need to reassess MT)

# Set the default AVX (vectorization) choice
ifeq ($(AVX),)
  ifeq ($(UNAME_P),ppc64le)
    ###override AVX = none
    override AVX = sse4
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

# Set the default HRDCOD (hardcode cIPC and cIPD physics parameters?) choice
ifeq ($(HRDCOD),)
  override HRDCOD = 0
endif

# Set the default RNDGEN (random number generator) choice
ifeq ($(NVCC),)
  override RNDGEN = hasNoCurand
else ifeq ($(RNDGEN),)
  override RNDGEN = hasCurand
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
else
  ifeq ($(AVX),sse4)
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
CXXFLAGS+= $(AVXFLAGS)

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
$(info FPTYPE=$(FPTYPE))
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE
  CUFLAGS  += -DMGONGPU_FPTYPE_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT
  CUFLAGS  += -DMGONGPU_FPTYPE_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'f' and 'd' are supported)
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
  CXXFLAGS += -DMGONGPU_HARDCODE_CIPC
  CUFLAGS  += -DMGONGPU_HARDCODE_CIPC
else ifneq ($(HRDCOD),0)
  $(error Unknown HRDCOD='$(HRDCOD)': only '0' and '1' are supported)
endif

# Set the build flags appropriate to each RNDGEN choice (example: "make RNDGEN=hasNoCurand")
$(info RNDGEN=$(RNDGEN))
ifeq ($(RNDGEN),hasNoCurand)
  CXXFLAGS += -DMGONGPU_HAS_NO_CURAND
else ifneq ($(RNDGEN),hasCurand)
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

cxx_main=$(BUILDDIR)/check.exe
cu_main=$(BUILDDIR)/gcheck.exe
testmain=$(BUILDDIR)/runTest.exe

# First target (default goal)
all.$(TAG): $(BUILDDIR)/.build.$(TAG) ../../src/$(BUILDDIR)/.build.$(TAG) $(cu_main) $(cxx_main) $(testmain)

# Target (and build options): debug
MAKEDEBUG=
debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: CUOPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

# Target: tag-specific build lockfiles
override oldtagsb=`if [ -d $(BUILDDIR) ]; then find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)' -exec echo $(shell pwd)/{} \; ; fi`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtagsb)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtagsb)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

../../src/$(BUILDDIR)/.build.$(TAG):
	$(MAKE) -C ../../src $(MAKEDEBUG)

#-------------------------------------------------------------------------------

# Generic target and build rules: objects from CUDA compilation
$(BUILDDIR)/%%.o : %%.cu *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) --compiler-options '-fPIC' -c $< -o $@

$(BUILDDIR)/%%_cu.o : %%.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) --compiler-options '-fPIC' -c -x cu $< -o $@

# Generic target and build rules: objects from C++ compilation
$(BUILDDIR)/%%.o : %%.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -fPIC -c $< -o $@

# Apply special build flags only to CrossSectionKernel.cc and gCrossSectionKernel.cu (no fast math, see #117)
$(BUILDDIR)/CrossSectionKernels.o: CXXFLAGS += -fno-fast-math
$(BUILDDIR)/gCrossSectionKernels.o: CUFLAGS += -Xcompiler -fno-fast-math

#### Apply special build flags only to CPPProcess.cc (-flto)
###$(BUILDDIR)/CPPProcess.o: CXXFLAGS += -flto

#### Apply special build flags only to CPPProcess.cc (AVXFLAGS)
###$(BUILDDIR)/CPPProcess.o: CXXFLAGS += $(AVXFLAGS)

#-------------------------------------------------------------------------------

# Target (and build rules): common (src) library
commonlib : $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so

$(LIBDIR)/lib$(MG5AMC_COMMONLIB).so: ../../src/*.h ../../src/*.cc
	$(MAKE) -C ../../src $(MAKEDEBUG)

#-------------------------------------------------------------------------------

MG5AMC_CXXLIB = mg5amc_cxx
cxx_objects=$(BUILDDIR)/CPPProcess.o $(BUILDDIR)/RandomNumberKernels.o $(BUILDDIR)/RamboSamplingKernels.o $(BUILDDIR)/MatrixElementKernels.o $(BUILDDIR)/BridgeKernels.o $(BUILDDIR)/CrossSectionKernels.o

MG5AMC_CULIB = mg5amc_cu
cu_objects=$(BUILDDIR)/gCPPProcess.o $(BUILDDIR)/gRandomNumberKernels.o $(BUILDDIR)/gRamboSamplingKernels.o $(BUILDDIR)/gMatrixElementKernels.o $(BUILDDIR)/gBridgeKernels.o $(BUILDDIR)/gCrossSectionKernels.o

# Target (and build rules): CUDA and C++ shared libraries
$(LIBDIR)/lib$(MG5AMC_CULIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cu_objects)
	$(NVCC) --shared -o $@ $(cu_objects) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

$(LIBDIR)/lib$(MG5AMC_CXXLIB).so: $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects)
	$(CXX) -shared -o $@ $(cxx_objects) -L$(LIBDIR) -l$(MG5AMC_COMMONLIB)

#-------------------------------------------------------------------------------

# Target (and build rules): CUDA and C++ standalone executables
ifneq ($(shell $(CXX) --version | grep ^Intel),)
$(cu_main): LIBFLAGS += -lintlc # compile with icpx and link with nvcc (undefined reference to `_intel_fast_memcpy')
endif

$(cu_main): $(BUILDDIR)/gcheck_sa.o $(LIBDIR)/lib$(MG5AMC_CULIB).so
	$(NVCC) -o $@ $(BUILDDIR)/gcheck_sa.o $(CUARCHFLAGS) $(LIBFLAGS) $(CULIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CULIB)

$(cxx_main): $(BUILDDIR)/check_sa.o $(LIBDIR)/lib$(MG5AMC_CXXLIB).so
	$(CXX) -o $@ $(BUILDDIR)/check_sa.o $(OMPFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS) -L$(LIBDIR) -l$(MG5AMC_CXXLIB)

#-------------------------------------------------------------------------------

# Target (and build rules): test objects and test executable
$(BUILDDIR)/testxxx.o: $(GTESTLIBS)
$(BUILDDIR)/testxxx.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(BUILDDIR)/testxxx.o: testxxx_cc_ref.txt
$(testmain): $(BUILDDIR)/testxxx.o
$(testmain): cxx_objects += $(BUILDDIR)/testxxx.o # Comment out this line to skip the C++ test of xxx functions

ifneq ($(NVCC),)
$(BUILDDIR)/testxxx_cu.o: $(GTESTLIBS)
$(BUILDDIR)/testxxx_cu.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(BUILDDIR)/testxxx_cu.o: testxxx_cc_ref.txt
$(testmain): $(BUILDDIR)/testxxx_cu.o
$(testmain): cu_objects += $(BUILDDIR)/testxxx_cu.o # Comment out this line to skip the CUDA test of xxx functions
endif

$(BUILDDIR)/testmisc.o: $(GTESTLIBS)
$(BUILDDIR)/testmisc.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): $(BUILDDIR)/testmisc.o
$(testmain): cxx_objects += $(BUILDDIR)/testmisc.o # Comment out this line to skip the C++ miscellaneous tests

ifneq ($(NVCC),)
$(BUILDDIR)/testmisc_cu.o: $(GTESTLIBS)
$(BUILDDIR)/testmisc_cu.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): $(BUILDDIR)/testmisc_cu.o
$(testmain): cu_objects += $(BUILDDIR)/testmisc_cu.o # Comment out this line to skip the CUDA miscellaneous tests
endif

$(BUILDDIR)/runTest.o: $(GTESTLIBS)
$(BUILDDIR)/runTest.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): $(BUILDDIR)/runTest.o
$(testmain): cxx_objects += $(BUILDDIR)/runTest.o

ifneq ($(NVCC),)
$(BUILDDIR)/runTest_cu.o: $(GTESTLIBS)
$(BUILDDIR)/runTest_cu.o: INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): $(BUILDDIR)/runTest_cu.o
$(testmain): cu_objects  += $(BUILDDIR)/runTest_cu.o
endif

$(testmain): $(GTESTLIBS)
$(testmain): INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main
ifneq ($(shell $(CXX) --version | grep ^clang),)
$(testmain): LIBFLAGS += -L$(patsubst %%bin/clang++,%%lib,$(shell which $(CXX) | tail -1))
endif

ifeq ($(NVCC),) # link only runTest.o
$(testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects) $(GTESTLIBS)
	$(CXX) -o $@ $(cxx_objects) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)
else # link both runTest.o and runTest_cu.o
$(testmain): $(LIBDIR)/lib$(MG5AMC_COMMONLIB).so $(cxx_objects) $(cu_objects) $(GTESTLIBS)
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) -ldl $(LIBFLAGS) $(CULIBFLAGS) -lcuda -lgomp
endif

$(GTESTLIBS):
	$(MAKE) -C $(TESTDIR)

#-------------------------------------------------------------------------------

# Target (and build rules): install libraries and headers (for use by MadEvent in Fortran)
INSTALL_HEADERS=CPPProcess.h MatrixElementKernels.h MemoryBuffers.h checkCuda.h MemoryAccessMomenta.h MemoryAccessHelpers.h MemoryAccessVectors.h
INSTALL_INC_DIR=../../include
INSTALL_OBJECTS=CPPProcess.o gCPPProcess.o MatrixElementKernels.o gMatrixElementKernels.o

install: all.$(TAG) $(INSTALL_INC_DIR) $(addprefix $(INSTALL_INC_DIR)/, $(INSTALL_HEADERS)) $(addprefix $(LIBDIR)/, $(INSTALL_OBJECTS))

$(INSTALL_INC_DIR) :
	$(MAKE) -C ../../src install

$(INSTALL_INC_DIR)/%%.h : %%.h
	@if [ ! -d $(INSTALL_INC_DIR) ]; then mkdir $(INSTALL_INC_DIR); fi
	cp $< $@

$(LIBDIR)/%%.o : $(BUILDDIR)/%%.o
	cp $< $@

#-------------------------------------------------------------------------------

# Target: build all targets in all AVX modes (each AVX mode in a separate build directory)
# Split the avxall target into five separate targets to allow parallel 'make -j avxall' builds
avxnone:
	@echo
	make USEBUILDDIR=1 AVX=none

avxsse4:
	@echo
	make USEBUILDDIR=1 AVX=sse4

avxavx2:
	@echo
	make USEBUILDDIR=1 AVX=avx2

avx512y:
	@echo
	make USEBUILDDIR=1 AVX=512y

avx512z:
	@echo
	make USEBUILDDIR=1 AVX=512z

avxall: avxnone avxsse4 avxavx2 avx512y avx512z

# Target: install all targets in all AVX modes (each AVX mode in a separate build directory)
# Do not make the avxall-install target dependent on avxall and do not implement -j (do 'make avxall -j; make avxall-install' instead)
avxall-install:
	@echo
	make USEBUILDDIR=1 AVX=none install
	@echo
	make USEBUILDDIR=1 AVX=sse4 install
	@echo
	make USEBUILDDIR=1 AVX=avx2 install
	@echo
	make USEBUILDDIR=1 AVX=512y install
	@echo
	make USEBUILDDIR=1 AVX=512z install

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

# Target: clean the builds as well as the googletest installation
distclean: cleanall
	make -C $(TOOLSDIR) clean
	make -C $(TESTDIR) clean

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

#-------------------------------------------------------------------------------

# Target: check (run the C++ test executable)
# [NB THIS IS WHAT IS USED IN THE GITHUB CI!]
check: runTest

# Target: runTest (run the C++ test executable runTest.exe)
runTest: $(BUILDDIR)/runTest.exe
	LD_LIBRARY_PATH="$(LIBDIR):$(LD_LIBRARY_PATH)" $(BUILDDIR)/runTest.exe

# Target: runCheck (run the CUDA standalone executable gcheck.exe with a small number of events)
runCheck: $(BUILDDIR)/check.exe
	LD_LIBRARY_PATH="$(LIBDIR):$(LD_LIBRARY_PATH)" $(BUILDDIR)/check.exe -p 2 32 2

# Target: runGcheck (run the CUDA standalone executable gcheck.exe with a small number of events)
runGcheck: $(BUILDDIR)/gcheck.exe
	LD_LIBRARY_PATH="$(LIBDIR):$(LD_LIBRARY_PATH)" $(BUILDDIR)/gcheck.exe -p 2 32 2

# Target: memcheck (run the CUDA standalone executable gcheck.exe with a small number of events through cuda-memcheck)
memcheck: $(BUILDDIR)/gcheck.exe
	LD_LIBRARY_PATH="$(LIBDIR):$(LD_LIBRARY_PATH)" $(CUDA_HOME)/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes $(BUILDDIR)/gcheck.exe -p 2 32 2

#-------------------------------------------------------------------------------
