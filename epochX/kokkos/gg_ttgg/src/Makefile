KOKKOSPATH_CUDA ?= /home/jchilders/git/kokkos/install_v100
KOKKOSPATH_OMP ?= /home/jchilders/git/kokkos/install_omp
KOKKOSPATH_INTEL ?= $(KOKKOS_HOME)
KOKKOSPATH_HIP ?= /home/jchilders/git/kokkos/install_mi50

MODELSTR = sm
MODELLIB = model_$(MODELSTR)

CXXFLAGS=-O3 -ffast-math --std=c++17

CUDA_ARCH_NUM ?= 70
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_CUDA)/include -arch=compute_$(CUDA_ARCH_NUM) --expt-extended-lambda --expt-relaxed-constexpr -use_fast_math --openmp -lineinfo

CXX=g++
ICPX=icpx
CLANG=clang++
HIPCC=hipcc

OPENMP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_OMP)/include --openmp

INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -fiopenmp -fopenmp-targets=spir64_gen -Wno-parentheses -Wno-openmp-mapping -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-unknown-attributes

HIP_ARCH_NUM ?= gfx906
HIP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_HIP)/include -fopenmp -fno-gpu-rdc --amdgpu-target=$(HIP_ARCH_NUM)

SRCS = Parameters_$(MODELSTR).cc read_slha.cc
# HelAmps_$(MODELSTR).cc

LIBDIR=../lib
cuda_lib=$(LIBDIR)/lib$(MODELLIB)_cuda.a
openmp_lib=$(LIBDIR)/lib$(MODELLIB)_openmp.a
intel_lib=$(LIBDIR)/lib$(MODELLIB)_intel.a
hip_lib=$(LIBDIR)/lib$(MODELLIB)_hip.a
openmp_objects=$(SRCS:.cc=.openmp.o)
cuda_objects=$(SRCS:.cc=.cuda.o)
intel_objects=$(SRCS:.cc=.intel.o)
hip_objects=$(SRCS:.cc=.hip.o)

all: cuda openmp intel hip

cuda: $(cuda_lib)

openmp: $(openmp_lib)

intel: $(intel_lib)

hip: $(hip_lib)

%.openmp.o : %.cc %.h
	$(CXX) $(OPENMP_CXXFLAGS) -c $< -o $@

%.cuda.o : %.cc %.h
	$(NVCC) $(CUDA_CXXFLAGS) -c $< -o $@

%.intel.o : %.cc %.h
	$(ICPX) $(INTEL_CXXFLAGS) -c $< -o $@

%.hip.o : %.cc %.h
	$(HIPCC) $(HIP_CXXFLAGS) -c $< -o $@

$(cuda_lib): $(cuda_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

$(openmp_lib): $(openmp_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

$(intel_lib): $(intel_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

$(hip_lib): $(hip_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

clean:
	rm -f $(cuda_objects) $(openmp_objects) $(intel_objects) $(hip_objects)
	rm -f $(cuda_lib) $(openmp_lib) $(intel_lib) $(hip_lib)
