KOKKOSPATH_CUDA ?= $(KOKKOS_HOME)
KOKKOSPATH_OMP ?= $(KOKKOS_HOME)
KOKKOSPATH_INTEL ?= $(KOKKOS_HOME)
KOKKOSPATH_HIP ?= $(KOKKOS_HOME)

processid_short=$(shell basename $(CURDIR) | awk -F_ '{print $$(NF-1)"_"$$NF}')
KOKKOS_BRIDGELIB=mg5amc_$(processid_short)_kokkos

MODELSTR = sm
MODELLIB = model_$(MODELSTR)
HELAMP_H = ../../src/HelAmps_$(MODELSTR).h
PAR_H = ../../src/Parameters_$(MODELSTR).h

INCDIR=../../src
CXXFLAGS=-O3 -ffast-math -I$(INCDIR) --std=c++17
LDFLAGS=-lstdc++fs

CUDA_ARCH_NUM ?= 70
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_CUDA)/include -arch=compute_$(CUDA_ARCH_NUM) --expt-extended-lambda --expt-relaxed-constexpr -use_fast_math --openmp -lineinfo
CUDA_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_CUDA)/lib64/libkokkoscore.a -lnvToolsExt --openmp

CXX ?= g++
ICPX ?= icpx
CLANG ?= clang++
HIPCC ?= hipcc

OPENMP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_OMP)/include --openmp
OPENMP_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_OMP)/lib64/libkokkoscore.a -ldl --openmp

# INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -I/soft/restricted/CNDA/sdk/2021.04.30.001/oneapi/compiler/latest/linux/include/sycl -fiopenmp -fopenmp-targets=spir64 -Wno-parentheses -Wno-openmp-mapping
# INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a -fiopenmp -fopenmp-targets=spir64 -L/soft/restricted/CNDA/sdk/2021.04.30.001/oneapi/compiler/latest/linux/lib/ -lsycl

INTEL_BACKEND ?= gen9
INTEL_SYCL ?= 0
ifeq ($(INTEL_SYCL),1)
# build for SYCL
INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -fiopenmp -fsycl -fsycl-targets=spir64_gen -Wno-parentheses -Wno-openmp-mapping -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-unknown-attributes
INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a $(KOKKOSPATH_INTEL)/lib64/libkokkoscontainers.a -fiopenmp -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend "-device ${INTEL_BACKEND}"
else
# build for OPENMP
INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -fiopenmp -fopenmp-targets=spir64_gen -Wno-parentheses -Wno-openmp-mapping -Wno-deprecated-declarations -Wno-tautological-constant-compare -Wno-unknown-attributes
INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a $(KOKKOSPATH_INTEL)/lib64/libkokkoscontainers.a -fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device ${INTEL_BACKEND}"
endif

HIP_ARCH_NUM ?= gfx906
HIP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_HIP)/include -fopenmp -fno-gpu-rdc --amdgpu-target=$(HIP_ARCH_NUM)
HIP_LDFLAGS=$(LDFLAGS) -L $(KOKKOSPATH_HIP)/lib64  -lkokkoscore -ldl -fopenmp 

FILES = check_sa.ext
CHECK_H = CPPProcess.h $(HELAMP_H) $(PAR_H) ../../src/rambo.h ../../src/random_generator.h CalcMean.h

cuda_exe=ccheck.exe
cuda_objects=$(FILES:.ext=.cuda.o)
cuda_libmodel=../../lib/lib$(MODELLIB)_cuda.a
openmp_exe=ocheck.exe
openmp_objects=$(FILES:.ext=.openmp.o)
openmp_libmodel=../../lib/lib$(MODELLIB)_openmp.a
intel_exe=icheck.exe
intel_objects=$(FILES:.ext=.intel.o)
intel_libmodel=../../lib/lib$(MODELLIB)_intel.a
hip_exe=hcheck.exe
hip_objects=$(FILES:.ext=.hip.o)
hip_libmodel=../../lib/lib$(MODELLIB)_hip.a
 
BRIDGE_SRCS = fbridge.ext
cuda_bridge_objects=$(BRIDGE_SRCS:.ext=.cuda.o)
openmp_bridge_objects=$(BRIDGE_SRCS:.ext=.openmp.o)
intel_bridge_objects=$(BRIDGE_SRCS:.ext=.intel.o)
hip_bridge_objects=$(BRIDGE_SRCS:.ext=.hip.o)

LIBDIR=../../lib
cuda_bridge_lib=lib$(MODELLIB)_bridge_cuda.a
openmp_bridge_lib=lib$(MODELLIB)_bridge_openmp.a
intel_bridge_lib=lib$(MODELLIB)_bridge_intel.a
hip_bridge_lib=lib$(MODELLIB)_bridge_hip.a


all: cuda openmp intel hip

cuda:
	make -C ../../src -f kokkos_src.mk cuda
	make $(cuda_exe)
openmp: 
	make -C ../../src -f kokkos_src.mk openmp
	make $(openmp_exe)
intel:
	INTEL_SYCL=$(INTEL_SYCL) make -C ../../src -f kokkos_src.mk intel
	make $(intel_exe)
hip:
	make -C ../../src -f kokkos_src.mk hip
	make $(hip_exe)

# compile object files

%.openmp.o : %.cc $(CHECK_H)
	$(CXX) $(OPENMP_CXXFLAGS) -c $< -o $@

%.cuda.o : %.cc $(CHECK_H)
	$(NVCC) $(CUDA_CXXFLAGS) -c $< -o $@

%.intel.o : %.cc $(CHECK_H)
	$(ICPX) $(INTEL_CXXFLAGS) -c $< -o $@

%.hip.o : %.cc $(CHECK_H)
	$(HIPCC) $(HIP_CXXFLAGS) -c $< -o $@

$(cuda_exe): $(cuda_objects) $(cuda_libmodel)
	$(NVCC) $(cuda_objects) -o $@ $(CUDA_LDFLAGS) $(cuda_libmodel)

$(openmp_exe): $(openmp_objects) $(openmp_libmodel)
	$(CXX) $(openmp_objects) -o $@ $(openmp_libmodel) $(OPENMP_LDFLAGS)

$(intel_exe): $(intel_objects) $(intel_libmodel)
	$(ICPX) $(intel_objects) -o $@ $(intel_libmodel) $(INTEL_LDFLAGS)

$(hip_exe): $(hip_objects) $(hip_libmodel)
	$(HIPCC) $(hip_objects) -o $@ $(HIP_LDFLAGS) $(hip_libmodel)


$(LIBDIR)/$(cuda_bridge_lib): $(cuda_bridge_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@
	ln -s $(cuda_bridge_lib) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a

$(LIBDIR)/$(openmp_bridge_lib): $(openmp_bridge_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@
	ln -s $(openmp_bridge_lib) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a

$(LIBDIR)/$(intel_bridge_lib): $(intel_bridge_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@
	ln -s $(intel_bridge_lib) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a

$(LIBDIR)/$(hip_bridge_lib): $(hip_bridge_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@
	ln -s $(hip_bridge_lib) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a

clean:
	make -C ../../src -f kokkos_src.mk clean
	rm -f $(cuda_objects) $(cuda_exe) $(openmp_objects) $(openmp_exe) $(intel_exe) $(intel_objects) $(hip_exe) $(hip_objects)
	rm -f $(cuda_bridge_objects) $(openmp_bridge_objects) $(intel_bridge_objects) $(hip_bridge_objects)
	rm -f $(LIBDIR)/$(cuda_bridge_lib) $(LIBDIR)/$(openmp_bridge_lib) $(LIBDIR)/$(intel_bridge_lib) $(LIBDIR)/$(hip_bridge_lib) $(LIBDIR)/lib$(KOKKOS_BRIDGELIB).a
