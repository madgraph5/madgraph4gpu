LIBDIR   = ../../lib/build.$(TAG)
TOOLSDIR = ../../../../../tools
TESTDIR  = ../../../../../test
INCFLAGS = -I. -I../../src -I$(TOOLSDIR)
MODELLIB = model_sm
OPTFLAGS = -O3 # this ends up in CUFLAGS too (should it?), cannot add -Ofast or -ffast-math here
OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

# AVX choice (example: "make AVX=none")
AVX ?= avx2
ifeq ($(AVX),avx2)
CXXFLAGS+= -march=core-avx2 # AVX2 (DEFAULT!)
else ifeq ($(AVX),avx512)
CXXFLAGS+= -march=native # AVX512F
else ifneq ($(AVX),none)
$(error Unknown AVX='$(AVX)': only 'avx2', 'avx512' and 'none' are supported)
endif

# Build tag (path to the build directory)
TAG = $(AVX)

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
  #CUARCHNUM=61 #(For Pascal Architecture Cards)
  USE_NVTX?=-DUSE_NVTX
  CUARCHFLAGS= -arch=compute_$(CUARCHNUM)
  # CUARCHFLAGS= -gencode arch=compute_$(CUARCHNUM),code=sm_$(CUARCHNUM)
  CUINC       = -I$(CUDA_HOME)/include
  CULIBFLAGS  = -L$(CUDA_HOME)/lib64 -lcuda -lcurand
  CUOPTFLAGS  = -lineinfo
  CUFLAGS= $(OPTFLAGS) $(CUOPTFLAGS) -std=c++14 $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math $(MGONGPU_CONFIG)
  # Without -maxrregcount: baseline throughput: 6.5E8 (16384 32 12) up to 7.3E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 160 # improves throughput: 6.9E8 (16384 32 12) up to 7.7E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 128 # improves throughput: 7.3E8 (16384 32 12) up to 7.6E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 96 # degrades throughput: 4.1E8 (16384 32 12) up to 4.5E8 (65536 128 12)
  ###CUFLAGS+= --maxrregcount 64 # degrades throughput: 1.7E8 (16384 32 12) flat at 1.7E8 (65536 128 12)
  cu_main     = build.$(TAG)/gcheck.exe
  cu_objects  = build.$(TAG)/gCPPProcess.o
else
  # No cuda. Switch cuda compilation off and go to common random numbers in C++
  NVCC       := $(warning CUDA_HOME is not set or is invalid. Export CUDA_HOME to compile with cuda)
  USE_NVTX   :=
  CULIBFLAGS :=
  ifndef MGONGPU_CONFIG
    export MGONGPU_CONFIG = -DMGONGPU_COMMONRAND_ONHOST
  endif
endif

GTESTLIBDIR = $(TESTDIR)/googletest/build/lib/
GTESTLIBS   = $(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a

MAKEDEBUG=

cxx_main=build.$(TAG)/check.exe
cxx_objects=build.$(TAG)/CPPProcess.o

testmain=build.$(TAG)/runTest.exe

# Assuming uname is available, detect if architecture is power
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
    CUFLAGS+= -Xcompiler -mno-float128
endif

all.$(TAG): ../../src $(cu_main) $(cxx_main) $(testmain)

debug: OPTFLAGS   = -g -O0 -DDEBUG2
debug: CUOPTFLAGS = -G
debug: MAKEDEBUG := debug
debug: all.$(TAG)

$(LIBDIR)/lib$(MODELLIB).a:
	$(MAKE) -C ../../src AVX=$(AVX) $(MAKEDEBUG)

build.$(TAG)/gcheck.o: gcheck.cu *.h ../../src/*.h ../../src/*.cu
	if [ ! -d build.$(TAG) ]; then mkdir build.$(TAG); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

build.$(TAG)/%.o : %.cu *.h ../../src/*.h
	if [ ! -d build.$(TAG) ]; then mkdir build.$(TAG); fi
	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

build.$(TAG)/%.o : %.cc *.h ../../src/*.h
	if [ ! -d build.$(TAG) ]; then mkdir build.$(TAG); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -c $< -o $@

$(cu_main): build.$(TAG)/gcheck.o $(LIBDIR)/lib$(MODELLIB).a $(cu_objects)
	$(NVCC) $< -o $@ $(cu_objects) $(CUARCHFLAGS) $(LIBFLAGS) $(CULIBFLAGS)

$(cxx_main): build.$(TAG)/check.o $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects)
	$(CXX) $< -o $@ $(cxx_objects) $(CPPFLAGS) $(CXXFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)

$(testmain): INCFLAGS += -I$(TESTDIR)/googletest/googletest/include
$(testmain): INCFLAGS += -I$(TESTDIR)/include
$(testmain): LIBFLAGS += -L$(GTESTLIBDIR) -lgtest -lgtest_main
ifeq ($(NVCC),)
$(testmain): $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(GTESTLIBS) build.$(TAG)/runTest.o $(TESTDIR)/src/MadgraphTest.o $(TESTDIR)/include/*.h
	$(CXX) -o $@ $(cxx_objects) build.$(TAG)/runTest.o $(TESTDIR)/src/MadgraphTest.o $(CPPFLAGS) $(CXXFLAGS) -ldl -pthread $(LIBFLAGS) $(CULIBFLAGS)
else
$(testmain) build.$(TAG)/runTest_cu.o &: runTest.cc $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(cu_objects) $(GTESTLIBS) $(TESTDIR)/include/*.h
	$(NVCC) -o build.$(TAG)/runTest_cu.o -c -x cu runTest.cc $(CPPFLAGS) $(CUFLAGS)
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) build.$(TAG)/runTest_cu.o $(CPPFLAGS) $(CUFLAGS) -ldl $(LIBFLAGS) $(CULIBFLAGS) -lcuda -lgomp
endif

$(GTESTLIBS):
	$(MAKE) -C $(TESTDIR)

check: $(testmain)
	$(testmain)

.PHONY: clean

clean:
	make -C ../../src AVX=$(AVX) clean
	rm -rf build.$(TAG)
#	rm -f build.$(TAG)/*.o build.$(TAG)/*.exe

distclean: clean
	make -C $(TOOLSDIR) clean
	make -C $(TESTDIR) clean

memcheck: $(cu_main)
	/usr/local/cuda/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes build.$(TAG)/gcheck.exe 2 32 2

perf: force
	make clean && make
	time build.$(TAG)/gcheck.exe -p 16348 32 12 && date

test: force
	build.$(TAG)/gcheck.exe -v 1 32 1

info:
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
