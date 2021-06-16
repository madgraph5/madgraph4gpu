INCFLAGS = -I.
OPTFLAGS = -O3
OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
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

# Enable ccache if USECCACHE=1
ifeq ($(USECCACHE)$(shell echo $(CXX) | grep ccache),1)
override CXX:=ccache $(CXX)
endif
ifeq ($(USECCACHE)$(shell echo $(AR) | grep ccache),1)
override AR:=ccache $(AR)
endif
ifeq ($(USECCACHE)$(shell echo $(NVCC) | grep ccache),1)
override NVCC:=ccache $(NVCC)
endif

# Assuming uname is available, detect if architecture is power
#UNAME_P := $(shell uname -p)
#ifeq ($(UNAME_P),ppc64le)
#  CUFLAGS+= -Xcompiler -mno-float128
#endif

target=$(LIBDIR)/libmodel_sm.a
cxx_objects=$(addprefix $(BUILDDIR)/, Parameters_sm.o read_slha.o rambo.o)
cu_objects= # NB grambo.cu must be included by gcheck.cu (no rdc)

all.$(TAG): $(BUILDDIR)/.build.tag_$(TAG) $(target)

override oldtags=`find $(BUILDDIR) -maxdepth 1 -name '.build.tag_*' ! -name '.build.tag_$(TAG)'`
$(BUILDDIR)/.build.tag_$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir $(BUILDDIR)"; mkdir $(BUILDDIR); fi
	@if [ "$(oldtags)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtags)\nPlease run 'make clean' first or consider using 'make USEBUILDDIR=1 AVX=$(AVX)'"; exit 1; fi
	@touch $(BUILDDIR)/.build.tag_$(TAG)

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
	rm -rf $(LIBDIR)
ifneq ($(BUILDDIR),.)
	rm -rf $(BUILDDIR)
else
	rm -f $(BUILDDIR)/.build.tag* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

cleanall:
	@echo
	make clean
	@echo
	rm -rf build.*
