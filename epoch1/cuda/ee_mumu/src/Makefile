INCFLAGS = -I.
OPTFLAGS = -O3
OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

# AVX choice (example: "make AVX=none")
ifneq ($(findstring clang++,$(CXX)),)
AVX := none
$(warning Using AVX='$(AVX)' for clang builds)
else ifneq ($(AVX),)
###$(info Using AVX='$(AVX)' according to user input)
else ifneq ($(shell grep -c avx512f /proc/cpuinfo),0)
AVX := avx512
###$(info Using AVX='$(AVX)' as no user input exists)
else
AVX := avx2
$(warning Using AVX='$(AVX)' as no user input exists and host does not support avx512f)
endif
###$(info AVX=$(AVX))

# Set the build flags appropriate to each AVX
ifeq ($(AVX),sse42)
override AVXFLAGS = -march=nehalem # SSE4.2
else ifeq ($(AVX),avx2)
override AVXFLAGS = -march=haswell # AVX2 (was -march=core-avx2 in previous versions of gcc)
else ifeq ($(AVX),avx512)
# NB: -mprefer-vector-width=256 is not needed (fptype_v assumes 256 anyway) but is a safer option
override AVXFLAGS = -march=skylake-avx512 -mprefer-vector-width=256 # AVX512F with 256 width (DEFAULT!)
else ifneq ($(AVX),none)
$(error Unknown AVX='$(AVX)': only 'sse42', 'avx2', 'avx512' and 'none' are supported)
endif

# For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations
CXXFLAGS+= $(AVXFLAGS)

# Build tag (defines target and path to the optional build directory)
override TAG = $(AVX)

# Build directory: current directory by default, or build.$(TAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
override BUILDDIR = build.$(TAG)
override LIBDIR   = ../lib/$(BUILDDIR)
else
override BUILDDIR = .
override LIBDIR   = ../lib
endif
###$(info BUILDDIR=$(BUILDDIR))
###$(info Building in BUILDDIR=$(BUILDDIR) for tag=$(TAG))

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
  CUINC = -I$(CUDA_HOME)/include/
endif

# Assuming uname is available, detect if architecture is power
#UNAME_P := $(shell uname -p)
#ifeq ($(UNAME_P),ppc64le)
#  CUFLAGS+= -Xcompiler -mno-float128
#endif

target=$(LIBDIR)/libmodel_sm.a
cxx_objects=$(addprefix $(BUILDDIR)/, Parameters_sm.o read_slha.o rambo.o)
cu_objects= # NB grambo.cu must be included by gcheck.cu (no rdc)

all.$(TAG): checktag $(target)

checktag:
	if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	@echo $(TAG) > $(BUILDDIR)/build.tag_new
	@if [ ! -e $(BUILDDIR)/build.tag ]; then mv -f $(BUILDDIR)/build.tag_new $(BUILDDIR)/build.tag; elif ! cmp -s $(BUILDDIR)/build.tag_new $(BUILDDIR)/build.tag; then echo -e "Cannot build for tag=$(TAG) as old build exists for tag=`cat $(BUILDDIR)/build.tag`\nPlease run 'make clean' first"; exit 1; fi

debug: OPTFLAGS = -g -O0 -DDEBUG2
#debug: CUFLAGS := $(filter-out -lineinfo,$(CUFLAGS))
#debug: CUFLAGS += -G
debug: $(target)

# NB: cuda includes are needed in the C++ code for curand.h
$(BUILDDIR)/%.o : %.cc *.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CUINC) -c $< -o $@

#$(BUILDDIR)/%.o : %.cu *.h
#	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
#	$(NVCC) $(CPPFLAGS) $(CUFLAGS) -c $< -o $@

$(target): $(cxx_objects) #$(cu_objects)
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) #$(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(BUILDDIR)/build.tag $(BUILDDIR)/build.tag_new
ifneq ($(BUILDDIR),.)
	rm -rf $(LIBDIR) $(BUILDDIR)
else
	rm -rf $(LIBDIR) $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif
