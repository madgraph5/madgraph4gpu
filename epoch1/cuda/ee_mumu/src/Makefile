INCFLAGS = -I.
OPTFLAGS = -O3
OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
LIBDIR   = ../lib/build.$(TAG)
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

# Export CUDA_HOME to select a cuda installation
ifndef CUDA_HOME
  NVCC ?= $(shell which nvcc 2>/dev/null)
  ifneq ($(NVCC),)
    # NVCC is in the PATH or set explicitly
    CUDA_HOME  = $(patsubst %bin/nvcc,%,$(NVCC))
    CUDA_HOME := $(warning No CUDA_HOME exported. Using "$(CUDA_HOME)") $(CUDA_HOME)
  endif
endif

ifdef CUDA_HOME
  CUINC = -I$(CUDA_HOME)/include
endif

# Assuming uname is available, detect if architecture is power
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
    CUFLAGS+= -Xcompiler -mno-float128
endif

target=$(LIBDIR)/libmodel_sm.a
cxx_objects=$(addprefix build.$(TAG)/, Parameters_sm.o read_slha.o rambo.o)
cu_objects= # NB grambo.cu must be included by gcheck.cu (no rdc)

all.$(TAG): $(target)

debug: OPTFLAGS = -g -O0 -DDEBUG2
debug: CUFLAGS := $(filter-out -lineinfo,$(CUFLAGS))
debug: CUFLAGS += -G
debug: $(target)

# NB: cuda includes are needed in the C++ code for curand.h
build.$(TAG)/%.o : %.cc *.h
	if [ ! -d build.$(TAG) ]; then mkdir build.$(TAG); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -c $< -o $@

#build.$(TAG)/%.o : %.cu *.h
#	if [ ! -d build.$(TAG) ]; then mkdir build.$(TAG); fi
#	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

$(target): $(cxx_objects) $(cu_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -rf build.$(TAG) $(LIBDIR)
#	rm -f $(cxx_objects) $(cu_objects) $(target)
