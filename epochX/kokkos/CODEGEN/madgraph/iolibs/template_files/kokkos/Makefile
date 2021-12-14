KOKKOSPATH_CUDA ?= /home/jchilders/git/kokkos/install_v100
KOKKOSPATH_OMP ?= /home/jchilders/git/kokkos/install_omp
KOKKOSPATH_INTEL ?= $(KOKKOS_HOME)
KOKKOSPATH_HIP ?= /home/jchilders/git/kokkos/install_mi50

MODELSTR = %(model)s
MODELLIB = model_$(MODELSTR)
HELAMP_H = ../../src/HelAmps_$(MODELSTR).h
PAR_H = ../../src/Parameters_$(MODELSTR).h

INCDIR=../../src
CXXFLAGS=-O3 -ffast-math -I$(INCDIR) --std=c++17 -lineinfo
LDFLAGS=

CUDA_ARCH_NUM ?= 70
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_CUDA)/include -arch=compute_$(CUDA_ARCH_NUM) --expt-extended-lambda --expt-relaxed-constexpr -use_fast_math --openmp
CUDA_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_CUDA)/lib64/libkokkoscore.a -lnvToolsExt --openmp

CXX=g++
IPCX=icpx
CLANG=clang++
HIPCC=hipcc

OPENMP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_OMP)/include --openmp
OPENMP_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_OMP)/lib64/libkokkoscore.a -ldl --openmp

# INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -I/soft/restricted/CNDA/sdk/2021.04.30.001/oneapi/compiler/latest/linux/include/sycl -fiopenmp -fopenmp-targets=spir64 -Wno-parentheses -Wno-openmp-mapping
# INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a -fiopenmp -fopenmp-targets=spir64 -L/soft/restricted/CNDA/sdk/2021.04.30.001/oneapi/compiler/latest/linux/lib/ -lsycl

INTEL_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_INTEL)/include -fiopenmp -fopenmp-targets=spir64 -Wno-parentheses -Wno-openmp-mapping
INTEL_LDFLAGS=$(LDFLAGS)  $(KOKKOSPATH_INTEL)/lib64/libkokkoscore.a -fiopenmp -fopenmp-targets=spir64 

HIP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_HIP)/include -fopenmp -fno-gpu-rdc --amdgpu-target=gfx906
HIP_LDFLAGS=$(LDFLAGS) -L $(KOKKOSPATH_HIP)/lib64  -lkokkoscore -ldl -fopenmp 

FILES = check.ext
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


all: cuda openmp intel hip

cuda:
	make -C ../../src cuda
	make $(cuda_exe)
openmp: 
	make -C ../../src openmp
	make $(openmp_exe)
intel:
	make -C ../../src intel
	make $(intel_exe)
hip:
	make -C ../../src hip
	make $(hip_exe)

# compile object files

%%.openmp.o : %%.cpp $(CHECK_H)
	$(CXX) $(OPENMP_CXXFLAGS) -c $< -o $@

%%.cuda.o : %%.cpp $(CHECK_H)
	$(NVCC) $(CUDA_CXXFLAGS) -c $< -o $@

%%.intel.o : %%.cpp $(CHECK_H)
	$(ICPX) $(INTEL_CXXFLAGS) -c $< -o $@

%%.hip.o : %%.cpp $(CHECK_H)
	$(HIPCC) $(HIP_CXXFLAGS) -c $< -o $@

$(cuda_exe): $(cuda_objects) $(cuda_libmodel)
	$(NVCC) $(cuda_objects) -o $@ $(CUDA_LDFLAGS) $(cuda_libmodel)

$(openmp_exe): $(openmp_objects) $(openmp_libmodel)
	$(CXX) $(openmp_objects) -o $@ $(OPENMP_LDFLAGS) $(openmp_libmodel)

$(intel_exe): $(intel_objects) $(intel_libmodel)
	$(IPCX) $(intel_objects) -o $@ $(INTEL_LDFLAGS) $(intel_libmodel)

$(hip_exe): $(hip_objects) $(hip_libmodel)
	$(HIPCC) $(hip_objects) -o $@ $(HIP_LDFLAGS) $(hip_libmodel)


clean:
	make -C ../../src clean
	rm -f $(cuda_objects) $(cuda_exe) $(openmp_objects) $(openmp_exe) $(intel_exe) $(intel_objects) $(hip_exe) $(hip_objects)
