KOKKOSPATH_CUDA=/home/parton/git/kokkos/install_v100
KOKKOSPATH_OMP=/home/parton/git/kokkos/install_omp
KOKKOSPATH_IRIS=/home/psteinbr/shared/gold/install

INCDIR=../../src
CXXFLAGS=-O3 -ffast-math -I$(INCDIR) --std=c++17
LDFLAGS=

CUDA_ARCH_NUM=70
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_CUDA)/include -arch=compute_$(CUDA_ARCH_NUM) --expt-extended-lambda --expt-relaxed-constexpr -use_fast_math --openmp
CUDA_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_CUDA)/lib64/libkokkoscore.a -lnvToolsExt --openmp

CXX=g++
IPCX=icpx

OPENMP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_OMP)/include --openmp
OPENMP_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_OMP)/lib64/libkokkoscore.a -ldl --openmp

IRIS_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_IRIS)/include -fiopenmp -fopenmp-targets=spir64 
IRIS_LDFLAGS=$(LDFLAGS) -L $(KOKKOSPATH_IRIS)/lib64 -lkokkoscore -fiopenmp -fopenmp-targets=spir64

cuda_exe=ccheck.exe
cuda_obj=ccheck_sa.o
openmp_exe=ocheck.exe
openmp_obj=ocheck_sa.o
iris_exe=icheck.exe
iris_obj=icheck.o
libmodel_sa=../../lib/libmodel_sa.a
iris_libmodel_sa=../../lib/libmodel_sa_iris.a

all: cuda openmp iris


cuda:
	make -C ../../src
	make $(cuda_exe)
openmp: 
	make -C ../../src
	make $(openmp_exe)
iris:
	make -C ../../src iris
	make $(iris_exe)


$(cuda_obj): check_sa.cpp ../../src/rambo.h ../../src/random_generator.h CPPProcess.h
	$(NVCC) -c $< -o $@ $(CUDA_CXXFLAGS)

$(cuda_exe): $(cuda_obj) $(libmodel_sa)
	$(NVCC)  $(cuda_obj) -o $@ $(CUDA_LDFLAGS) $(libmodel_sa)

$(openmp_obj): check_sa.cpp ../../src/rambo.h ../../src/random_generator.h CPPProcess.h
	$(CXX) -c $< -o $@ $(OPENMP_CXXFLAGS)

$(openmp_exe): $(openmp_obj) $(libmodel_sa)
	$(CXX) $(openmp_obj) -o $@ $(OPENMP_LDFLAGS) $(libmodel_sa)

$(iris_obj): check_sa.cpp ../../src/rambo.h ../../src/random_generator.h CPPProcess.h
	$(IPCX) -c $< -o $@ $(IRIS_CXXFLAGS)

$(iris_exe): $(iris_obj) $(iris_libmodel_sa)
	$(IPCX) $(iris_obj) -o $@ $(IRIS_LDFLAGS) $(iris_libmodel_sa)


.PHONY: clean

clean:
	make -C ../../src clean
	rm -f $(cuda_obj) $(cuda_exe) $(openmp_obj) $(openmp_exe) $(iris_exe) $(iris_obj)

memcheck: check
	/usr/local/cuda/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes ./check.exe 5 5 10

perf: force
	make clean && make
	time ./check.exe -p 1024 16 384 && date

test: force
	./check.exe -v 1 1 10

force:

#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
