TOOLSDIR = ../../../../../tools
TESTDIR  = ../../../../../test
INCFLAGS = -I. -I../../src -I$(TOOLSDIR)
MODELLIB = model_sm
OPTFLAGS = -O3 # this ends up in CUFLAGS too (should it?), cannot add -Ofast or -ffast-math here
ifneq ($(shell $(CXX) --version | grep ^Intel),)
OMPFLAGS?=
else
OMPFLAGS?= -fopenmp
endif
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
###CXXFLAGS+= -g # FOR DEBUGGING ONLY
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

# Assuming uname is available, detect if architecture is PowerPC
UNAME_P := $(shell uname -p)

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

# Set the default RNDGEN (random number generator) choice
ifeq ($(RNDGEN),)
  override RNDGEN = common
endif

# If CUDA_HOME is not set, try to set it from the location of nvcc
ifndef CUDA_HOME
  CUDA_HOME = $(patsubst %bin/nvcc,%,$(shell which nvcc 2>/dev/null))
  $(warning CUDA_HOME was not set: using "$(CUDA_HOME)")
endif

# Set NVCC as $(CUDA_HOME)/bin/nvcc if it exists
ifneq ($(wildcard $(CUDA_HOME)/bin/nvcc),)
  NVCC = $(CUDA_HOME)/bin/nvcc
  CUARCHNUM=70
  ###CUARCHNUM=61 # (For Pascal Architecture Cards)
  USE_NVTX ?=-DUSE_NVTX
  CUARCHFLAGS = -arch=compute_$(CUARCHNUM)
  ###CUARCHFLAGS = -gencode arch=compute_$(CUARCHNUM),code=sm_$(CUARCHNUM)
  CUINC       = -I$(CUDA_HOME)/include/
  CULIBFLAGS  = -L$(CUDA_HOME)/lib64/ -lcurand # NB: -lcuda is not needed here!
  CUOPTFLAGS  = -lineinfo
  CUFLAGS     = $(OPTFLAGS) $(CUOPTFLAGS) -std=c++14 $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)
  cu_main     = $(BUILDDIR)/gcheck.exe
  cu_objects  = $(BUILDDIR)/gCPPProcess.o

  # See if we can build a CUPLA version targetting cuda
  ifneq ($(wildcard $(CUPLA_ROOT)/include/cuda_to_cupla.hpp),)
    ifeq ($(wildcard $(ALPAKA_ROOT)/include/alpaka/alpaka.hpp),)
      ALPAKA_ROOT := $(CUPLA_ROOT)/alpaka
    endif
    alp_main = $(BUILDDIR)/alpcheck.exe
    alp_objects = $(BUILDDIR)/alpCPPProcess.o
    cupla_objects = $(BUILDDIR)/cupla/common.o $(BUILDDIR)/cupla/device.o  $(BUILDDIR)/cupla/event.o $(BUILDDIR)/cupla/memory.o $(BUILDDIR)/cupla/stream.o $(BUILDDIR)/cupla/manager/Driver.o
    ALPINC := -I$(ALPAKA_ROOT)/include -I$(BOOSTINC)
    CUPLAINC := -I$(CUPLA_ROOT)/include
    CUPLASRC := $(CUPLA_ROOT)/src
    ALPFLAGS := -DALPAKA -DALPAKA_ACC_GPU_CUDA_ENABLED --expt-relaxed-constexpr
  endif

else
  # No cuda. Switch cuda compilation off and go to common random numbers in C++
  $(warning CUDA_HOME is not set or is invalid. Export CUDA_HOME to compile with cuda)
  NVCC       := 
  USE_NVTX   :=
  CULIBFLAGS :=
  override RNDGEN = common
endif

# Set the host compiler for nvcc
CUFLAGS += -ccbin $(shell which $(subst ccache ,,$(CXX)))

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

# Set the build flags appropriate to each RNDGEN choice (example: "make RNDGEN=common")
$(info RNDGEN=$(RNDGEN))
ifeq ($(RNDGEN),curdev)
  CXXFLAGS += -DMGONGPU_CURAND_ONDEVICE
  CUFLAGS  += -DMGONGPU_CURAND_ONDEVICE
else ifeq ($(RNDGEN),curhst)
  CXXFLAGS += -DMGONGPU_CURAND_ONHOST
  CUFLAGS  += -DMGONGPU_CURAND_ONHOST
else ifeq ($(RNDGEN),common)
  CXXFLAGS += -DMGONGPU_COMMONRAND_ONHOST
  CUFLAGS  += -DMGONGPU_COMMONRAND_ONHOST
else
  $(error Unknown RNDGEN='$(RNDGEN)': only 'curdev', 'curhst' and 'common' are supported)
endif

# Export AVX, FPTYPE, HELINL, RNDGEN so that it is not necessary to pass them to the src Makefile too
export AVX
export FPTYPE
export HELINL
export RNDGEN

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
override DIRTAG = $(AVX)_$(FPTYPE)_inl$(HELINL)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of CUDA and no-CUDA environment builds with different random number generators)
override TAG = $(AVX)_$(FPTYPE)_inl$(HELINL)_$(RNDGEN)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
  override LIBDIR   = ../../lib/$(BUILDDIR)
else
  override BUILDDIR = .
  override LIBDIR   = ../../lib
endif
###$(info BUILDDIR=$(BUILDDIR))
$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

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

GTESTLIBDIR = $(TESTDIR)/googletest/build/lib/
GTESTLIBS   = $(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a

MAKEDEBUG=

cxx_main=$(BUILDDIR)/check.exe
cxx_objects=$(BUILDDIR)/CPPProcess.o

testmain=$(BUILDDIR)/runTest.exe

all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(LIBDIR)/lib$(MODELLIB).a $(cu_main) $(alp_main) $(cxx_main) $(testmain)

override oldtags=`find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)'`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	@if [ "$(oldtags)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtags)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: CUOPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

$(LIBDIR)/lib$(MODELLIB).a: ../../src/*.h ../../src/*.cc
	$(MAKE) -C ../../src $(MAKEDEBUG)

$(BUILDDIR)/gcheck_sa.o: gcheck_sa.cu *.h ../../src/*.h # ../../src/*.cu
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

$(BUILDDIR)/alpcheck_sa.o: gcheck_sa.cu *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(LIBFLAGS) $(CULIBFLAGS) $(CUFLAGS) $(ALPFLAGS) $(ALPINC) $(CUPLAINC) -dc $< -o $@

$(BUILDDIR)/CPPProcess.o : ../../src/HelAmps_sm.cc
$(BUILDDIR)/gCPPProcess.o : ../../src/HelAmps_sm.cc

$(BUILDDIR)/%.o : %.cu *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

# Apply special build flags only to CPPProcess.cc (-flto)
#$(BUILDDIR)/CPPProcess.o : CPPProcess.cc *.h ../../src/*.h
#	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -flto $(CUINC) -c $< -o $@

### Apply special build flags only to CPPProcess.cc (AVXFLAGS)
###$(BUILDDIR)/CPPProcess.o : CPPProcess.cc *.h ../../src/*.h
###	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
###	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(AVXFLAGS) $(CUINC) -c $< -o $@

$(BUILDDIR)/%.o : %.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir -p $(BUILDDIR)"; mkdir -p $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -c $< -o $@

$(BUILDDIR)/cupla/%.o: $(CUPLASRC)/%.cpp
	@if [ ! -d $(BUILDDIR)/cupla/manager ]; then mkdir -p $(BUILDDIR)/cupla/manager; fi
	$(NVCC) -x cu -dc $< -o $@ $(CULIBFLAGS) $(ALPFLAGS) $(ALPINC) $(CUPLAINC)

$(cu_main): $(BUILDDIR)/gcheck_sa.o $(LIBDIR)/lib$(MODELLIB).a $(cu_objects)
	$(NVCC) $< -o $@ $(cu_objects) $(CUARCHFLAGS) $(LIBFLAGS) $(CULIBFLAGS)

$(BUILDDIR)/alpCPPProcess.o: gCPPProcess.cu ../../src/HelAmps_sm.cc
	$(NVCC) $(CULIBFLAGS) $(LIBFLAGS) $(CUFLAGS) $(ALPFLAGS) $(ALPINC) $(CUPLAINC) -dc $< -o $@

$(alp_main): $(BUILDDIR)/alpcheck_sa.o $(LIBDIR)/lib$(MODELLIB).a $(alp_objects) $(cupla_objects)
	$(NVCC) $< -o $@ $(alp_objects) $(cupla_objects) $(CUFLAGS) $(LIBFLAGS) $(CULIBFLAGS)

$(cxx_main): $(BUILDDIR)/check_sa.o $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects)
	$(CXX) $< -o $@ $(cxx_objects) $(CPPFLAGS) $(CXXFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)

$(BUILDDIR)/runTest.o: $(GTESTLIBS)
$(testmain): $(GTESTLIBS)
$(testmain): INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main
$(testmain): $(BUILDDIR)/runTest.o $(BUILDDIR)/MadgraphTest.o
$(testmain): cxx_objects += $(BUILDDIR)/runTest.o $(BUILDDIR)/MadgraphTest.o
$(testmain): alp_objects += $(BUILDDIR)/runTest_alp.o
ifneq ($(shell $(CXX) --version | grep ^clang),)
$(testmain): LIBFLAGS += -L$(patsubst %bin/clang++,%lib,$(shell which $(CXX) | tail -1))
endif

$(BUILDDIR)/testxxx.o:   $(GTESTLIBS)
$(BUILDDIR)/testxxx.o:   testxxx_cc_ref.txt
$(testmain): $(BUILDDIR)/testxxx.o
$(testmain): cxx_objects += $(BUILDDIR)/testxxx.o # Comment out this line to skip the test of xxx functions

ifeq ($(NVCC),)
# Link only runTest.o
$(testmain): $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(GTESTLIBS) 
	$(CXX) -o $@ $(cxx_objects) $(CPPFLAGS) $(CXXFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)
else ifneq ($(wildcard $(CUPLA_ROOT)/include/cuda_to_cupla.hpp),)
# Link both runTest.o and runTest_alp.o
# (todo? avoid multiple targets and '&', this needs the latest make 4.3, see https://stackoverflow.com/a/60232515)
$(testmain) $(BUILDDIR)/runTest_alp.o &: runTest.cc $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(alp_objects) $(cupla_objects) $(GTESTLIBS)
	$(NVCC) -o $(BUILDDIR)/runTest_alp.o -dc -x cu runTest.cc $(CUFLAGS) $(ALPFLAGS) $(ALPINC) $(CUPLAINC)
	$(NVCC) -o $@ $(cxx_objects) $(alp_objects) $(cupla_objects) $(CPPFLAGS) $(CUFLAGS) -ldl $(LIBFLAGS) $(CULIBFLAGS) -lcuda -lgomp
else
# Link both runTest.o and runTest_cu.o
# (todo? avoid multiple targets and '&', this needs the latest make 4.3, see https://stackoverflow.com/a/60232515)
$(testmain) $(BUILDDIR)/runTest_cu.o &: runTest.cc $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(cu_objects) $(GTESTLIBS)
	$(NVCC) -o $(BUILDDIR)/runTest_cu.o -c -x cu runTest.cc $(CPPFLAGS) $(CUFLAGS)
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) $(CPPFLAGS) $(CUFLAGS) -ldl $(LIBFLAGS) $(CULIBFLAGS) -lcuda -lgomp
endif

$(GTESTLIBS):
	$(MAKE) -C $(TESTDIR)

check: $(testmain)
	$(testmain)

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

.PHONY: clean

clean:
	make -C ../../src clean
ifneq ($(BUILDDIR),.)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
	rm -rf cupla
endif

cleanall:
	@echo
	make clean
	@echo
	make -C ../../src cleanall
	rm -rf build.*	

distclean: cleanall
	make -C $(TOOLSDIR) clean
	make -C $(TESTDIR) clean

memcheck: $(cu_main)
	/usr/local/cuda/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes $(BUILDDIR)/gcheck.exe 2 32 2

perf: force
	make clean && make
	time $(BUILDDIR)/gcheck.exe -p 16348 32 12 && date

test: force
	$(BUILDDIR)/gcheck.exe -v 1 32 1

info:
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

force:

#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
