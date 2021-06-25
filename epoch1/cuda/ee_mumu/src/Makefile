INCFLAGS = -I.
OPTFLAGS = -O3
OMPFLAGS?= -fopenmp
CXXFLAGS = $(OPTFLAGS) -std=c++17 $(INCFLAGS) $(USE_NVTX) -Wall -Wshadow -Wextra $(OMPFLAGS) $(MGONGPU_CONFIG)
CXXFLAGS+= -ffast-math # see issue #117
###CXXFLAGS+= -Ofast # performance is not different from --fast-math
###CXXFLAGS+= -g # FOR DEBUGGING ONLY
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

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

# Set the host compiler for nvcc
override NVCC := $(NVCC) -ccbin $(CXX)

# Set the build flags appropriate to each AVX choice (example: "make AVX=none")
# [NB MGONGPU_PVW512 is needed because "-mprefer-vector-width=256" is not exposed in a macro]
# [See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96476]
###$(info AVX=$(AVX))
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
# For the moment, use AVXFLAGS everywhere: eventually, use them only in encapsulated implementations?
CXXFLAGS+= $(AVXFLAGS)

# Set the build flags appropriate to each FPTYPE choice (example: "make FPTYPE=f")
###$(info FPTYPE=$(FPTYPE))
ifeq ($(FPTYPE),d)
  CXXFLAGS += -DMGONGPU_FPTYPE_DOUBLE
  CUFLAGS  += -DMGONGPU_FPTYPE_DOUBLE
else ifeq ($(FPTYPE),f)
  CXXFLAGS += -DMGONGPU_FPTYPE_FLOAT
  CUFLAGS  += -DMGONGPU_FPTYPE_FLOAT
else
  $(error Unknown FPTYPE='$(FPTYPE)': only 'f' and 'd' are supported)
endif

# Set the build flags appropriate to each RNDGEN choice (example: "make RNDGEN=common")
###$(info RNDGEN=$(RNDGEN))
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

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
override DIRTAG = $(AVX)_$(FPTYPE)

# Build lockfile "full" tag (defines full specification of build options that cannot be intermixed)
# (Rationale: avoid mixing of CUDA and no-CUDA environment builds with different random number generators)
override TAG = $(AVX)_$(FPTYPE)_$(RNDGEN)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override BUILDDIR = build.$(DIRTAG)
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

all.$(TAG): $(BUILDDIR)/.build.$(TAG) $(target)

override oldtags=`find $(BUILDDIR) -maxdepth 1 -name '.build.*' ! -name '.build.$(TAG)'`
$(BUILDDIR)/.build.$(TAG):
	@if [ ! -d $(BUILDDIR) ]; then echo "mkdir $(BUILDDIR)"; mkdir $(BUILDDIR); fi
	@if [ "$(oldtags)" != "" ]; then echo -e "Cannot build for tag=$(TAG) as old builds exist for other tags:\n$(oldtags)\nPlease run 'make clean' first\nIf 'make clean' is not enough: run 'make clean USEBUILDDIR=1 AVX=$(AVX) FPTYPE=$(FPTYPE)' or 'make cleanall'"; exit 1; fi
	@touch $(BUILDDIR)/.build.$(TAG)

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
	rm -f $(BUILDDIR)/.build.* $(BUILDDIR)/*.o $(BUILDDIR)/*.exe
endif

cleanall:
	@echo
	make clean
	@echo
	rm -rf build.*
