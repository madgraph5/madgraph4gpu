KOKKOSPATH_CUDA=/home/jchilders/git/kokkos/install
KOKKOSPATH_OMP=/home/jchilders/git/kokkos/install_openmp

INCDIR=../../src
CXXFLAGS=-O3 -ffast-math -I$(INCDIR) --std=c++17 --openmp
LDFLAGS=--openmp

CUDA_ARCH_NUM=75
NVCC=$(KOKKOSPATH_CUDA)/bin/nvcc_wrapper
CUDA_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_CUDA)/include -arch=compute_$(CUDA_ARCH_NUM) --expt-extended-lambda --expt-relaxed-constexpr -use_fast_math
CUDA_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_CUDA)/lib/libkokkoscore.a -lnvToolsExt

CXX=g++

OPENMP_CXXFLAGS=$(CXXFLAGS) -I$(KOKKOSPATH_OMP)/include
OPENMP_LDFLAGS=$(LDFLAGS) $(KOKKOSPATH_OMP)/lib/libkokkoscore.a -ldl

cuda_exe=ccheck.exe
cuda_obj=ccheck_sa.o
openmp_exe=ocheck.exe
openmp_obj=ocheck_sa.o
libmodel_sa=../../lib/libmodel_sa.a

all: cuda openmp


cuda:
	make -C ../../src
	make $(cuda_exe)
openmp: 
	make -C ../../src
	make $(openmp_exe)


$(cuda_obj): check_sa.cpp ../../src/rambo.h ../../src/random_generator.h CPPProcess.h
	$(NVCC) -c $< -o $@ $(CUDA_CXXFLAGS)

$(cuda_exe): $(cuda_obj) $(libmodel_sa)
	$(NVCC)  $(cuda_obj) -o $@ $(CUDA_LDFLAGS) $(libmodel_sa)

$(openmp_obj): check_sa.cpp ../../src/rambo.h ../../src/random_generator.h CPPProcess.h
	$(CXX) -c $< -o $@ $(OPENMP_CXXFLAGS)

$(openmp_exe): $(openmp_obj) $(libmodel_sa)
	$(CXX) $(openmp_obj) -o $@ $(OPENMP_LDFLAGS) $(libmodel_sa)

.PHONY: clean

clean:
	make -C ../../src clean
	rm -f $(cuda_obj) $(cuda_exe) $(openmp_obj) $(openmp_exe)

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
