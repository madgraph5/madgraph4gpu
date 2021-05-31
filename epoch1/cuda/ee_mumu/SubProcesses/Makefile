TOOLSDIR = ../../../../../tools
TESTDIR  = ../../../../../test
INCFLAGS = -I. -I../../src -I$(TOOLSDIR)
MODELLIB = model_sm
OPTFLAGS = -O3 # this ends up in CUFLAGS too (should it?), cannot add -Ofast or -ffast-math here
#OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
CXXFLAGS+= -g # FOR DEBUGGING ONLY
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

# AVX choice (example: "make AVX=none")
ifneq ($(AVX),)
###$(info Using AVX='$(AVX)' according to user input)
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
###$(info AVX=$(AVX))

# Set the build flags appropriate to each AVX
# [NB MGONGPU_PVW512 is needed because "-mprefer-vector-width=256" is not exposed in a macro]
# [See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96476]
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

# For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations
CXXFLAGS+= $(AVXFLAGS)

# Build tag (defines target and path to the optional build directory)
override TAG = $(AVX)

# Build directory: current directory by default, or build.$(TAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
override BUILDDIR = build.$(TAG)
override LIBDIR   = ../../lib/$(BUILDDIR)
else
override BUILDDIR = .
override LIBDIR   = ../../lib
endif
###$(info BUILDDIR=$(BUILDDIR))
$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

# If CUDA_HOME is not set, try to set it from the location of nvcc
ifndef CUDA_HOME
  NVCC ?= $(shell which nvcc 2>/dev/null)
  ifneq ($(NVCC),)
    # NVCC is in the PATH or set explicitly
    CUDA_HOME  = $(patsubst %bin/nvcc,%,$(NVCC))
    CUDA_HOME := $(warning No CUDA_HOME exported. Using "$(CUDA_HOME)") $(CUDA_HOME)
  endif
endif

ifneq ($(wildcard $(CUDA_HOME)/bin/nvcc),)
  NVCC = $(CUDA_HOME)/bin/nvcc
  CUARCHNUM=70
  ###CUARCHNUM=61 # (For Pascal Architecture Cards)
  USE_NVTX ?=-DUSE_NVTX
  CUARCHFLAGS = -arch=compute_$(CUARCHNUM)
  ###CUARCHFLAGS = -gencode arch=compute_$(CUARCHNUM),code=sm_$(CUARCHNUM)
  CUINC       = -I$(CUDA_HOME)/include/
  CULIBFLAGS  = -L$(CUDA_HOME)/lib64/ -lcuda -lcurand
  CUOPTFLAGS  = -lineinfo
  CUFLAGS     = $(OPTFLAGS) $(CUOPTFLAGS) -std=c++14 $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math $(MGONGPU_CONFIG)
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)
  cu_main     = $(BUILDDIR)/gcheck.exe
  cu_objects  = $(BUILDDIR)/gCPPProcess.o
else
  # No cuda. Switch cuda compilation off and go to common random numbers in C++
  NVCC       := $(warning CUDA_HOME is not set or is invalid. Export CUDA_HOME to compile with cuda)
  USE_NVTX   :=
  CULIBFLAGS :=
  ifndef MGONGPU_CONFIG
    export MGONGPU_CONFIG = -DMGONGPU_COMMONRAND_ONHOST
  endif
endif

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
override CXX:=ccache $(CXX)
endif
ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
override AR:=ccache $(AR)
endif
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

# Assuming uname is available, detect if architecture is power
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
    CUFLAGS+= -Xcompiler -mno-float128
endif

all.$(TAG): ../../src/$(BUILDDIR)/.build.tag_$(TAG) $(BUILDDIR)/.build.tag_$(TAG) $(cu_main) $(cxx_main) $(testmain)

override oldtags=`find $(BUILDDIR) -maxdepth 1 -name '.build.tag_*' ! -name '.build.tag_$(TAG)'`
$(BUILDDIR)/.build.tag_$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir $(BUILDDIR)"; mkdir $(BUILDDIR); fi
	@if [ "$(oldtags)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtags)\nPlease run 'make clean' first or consider using 'make USEBUILDDIR=1 AVX=$(AVX)'"; exit 1; fi
	@touch $(BUILDDIR)/.build.tag_$(TAG)

debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: CUOPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

../../src/$(BUILDDIR)/.build.tag_$(TAG):
	$(MAKE) -C ../../src AVX=$(AVX) $(MAKEDEBUG)

$(LIBDIR)/lib$(MODELLIB).a: ../../src/*.h ../../src/*.cc
	$(MAKE) -C ../../src AVX=$(AVX) $(MAKEDEBUG)

$(BUILDDIR)/gcheck.o: gcheck.cu *.h ../../src/*.h ../../src/*.cu
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

$(BUILDDIR)/CPPProcess.o : ../../src/HelAmps_sm.cc
$(BUILDDIR)/gCPPProcess.o : ../../src/HelAmps_sm.cc

$(BUILDDIR)/%.o : %.cu *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

#$(BUILDDIR)/CPPProcess.o : CPPProcess.cc *.h ../../src/*.h
#	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(AVXFLAGS) $(CUINC) -c $< -o $@

$(BUILDDIR)/%.o : %.cc *.h ../../src/*.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -c $< -o $@

$(cu_main): $(BUILDDIR)/gcheck.o $(LIBDIR)/lib$(MODELLIB).a $(cu_objects)
	$(NVCC) $< -o $@ $(cu_objects) $(CUARCHFLAGS) $(LIBFLAGS) $(CULIBFLAGS)

$(cxx_main): $(BUILDDIR)/check.o $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects)
	$(CXX) $< -o $@ $(cxx_objects) $(CPPFLAGS) $(CXXFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)

$(BUILDDIR)/runTest.o: $(GTESTLIBS)
$(testmain): $(GTESTLIBS)
$(testmain): INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): INCFLAGS += -I$(TESTDIR)/include
$(testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main
$(testmain): $(BUILDDIR)/runTest.o $(TESTDIR)/src/MadgraphTest.o $(TESTDIR)/include/*.h
$(testmain): cxx_objects += $(BUILDDIR)/runTest.o $(TESTDIR)/src/MadgraphTest.o
$(testmain): cu_objects  += $(BUILDDIR)/runTest_cu.o
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
else
# Link both runTest.o and runTest_cu.o
# (todo? avoid multiple targets and '&', this needs the latest make 4.3, see https://stackoverflow.com/a/60232515)
$(testmain) $(BUILDDIR)/runTest_cu.o &: runTest.cc $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(cu_objects) $(GTESTLIBS)
	$(NVCC) -o $(BUILDDIR)/runTest_cu.o -c -x cu runTest.cc $(CPPFLAGS) $(CUFLAGS)
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) $(CPPFLAGS) $(CUFLAGS) -ldl $(LIBFLAGS) $(CULIBFLAGS) -lcuda -lgomp
endif

$(TESTDIR)/src/MadgraphTest.o: $(TESTDIR)/src/MadgraphTest.cc $(TESTDIR)/include/*.h
	$(MAKE) -C $(TESTDIR)

$(GTESTLIBS):
	$(MAKE) -C $(TESTDIR)

check: $(testmain)
	$(testmain)

.PHONY: clean

clean:
	make -C ../../src AVX=$(AVX) clean
	rm -f $(BUILDDIR)/.build.tag*
ifneq ($(BUILDDIR),.)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

avxall:
	@echo
	make USEBUILDDIR=1 AVX=none
	@echo
	make USEBUILDDIR=1 AVX=sse4
	@echo
	make USEBUILDDIR=1 AVX=avx2
	@echo
	make USEBUILDDIR=1 AVX=512y
	@echo
	make USEBUILDDIR=1 AVX=512z

cleanall:
	@echo
	make clean
	@echo
	make USEBUILDDIR=1 AVX=none clean; make -C ../../src USEBUILDDIR=1 AVX=none clean
	@echo
	make USEBUILDDIR=1 AVX=sse4 clean; make -C ../../src USEBUILDDIR=1 AVX=sse4 clean
	@echo
	make USEBUILDDIR=1 AVX=avx2 clean; make -C ../../src USEBUILDDIR=1 AVX=avx2 clean
	@echo
	make USEBUILDDIR=1 AVX=512y clean; make -C ../../src USEBUILDDIR=1 AVX=512y clean
	@echo
	make USEBUILDDIR=1 AVX=512z clean; make -C ../../src USEBUILDDIR=1 AVX=512z clean

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
ifeq ($(USECCACHE),1)
	ccache --version | head -1
	@echo ""
endif
ifneq ($(NVCC),)
	$(NVCC) --version
	@echo ""
endif
	$(CXX) --version

force:

#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
