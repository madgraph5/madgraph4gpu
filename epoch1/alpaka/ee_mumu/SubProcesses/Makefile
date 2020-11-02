LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -g -O3 -I. -I$(INCDIR) -std=c++17 -Wall -Wshadow -fopenmp
ALPINC=/data/dhsmith/alpaka-training/alpaka-install/include
CUPLAINC=/data/dhsmith/cupla/include

#CUPLAFLAGS=-DALPAKA_ACC_CPU_B_SEQ_T_THREADS_ENABLED
#CUPLAFLAGS=-DALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
#CUPLAFLAGS=-DALPAKA_ACC_CPU_B_OMP2_T_SEQ_ENABLED
CUPLAFLAGS=-DALPAKA_ACC_GPU_CUDA_ENABLED

CUPLASRC=/data/dhsmith/cupla/src
BOOSTINC=/usr/include/boost169
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB) -lpthread -fopenmp
CXX=g++
MAKEDEBUG=

cxx_main=check.exe

#cxx_objects=CPPProcess.o cupla/common.o cupla/device.o cupla/event.o cupla/memory.o cupla/stream.o cupla/manager/Driver.o
cxx_objects=cupla/common.o cupla/device.o cupla/event.o cupla/memory.o cupla/stream.o cupla/manager/Driver.o

all:
	cd ../../src && make
	mkdir -p cupla/manager
	make $(cxx_main)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG2
debug: MAKEDEBUG := debug
debug: $(cu_main)

$(LIBDIR)/lib$(MODELLIB).a:
	@cd ../../src && make $(MAKEDEBUG)

cupla/%.o: $(CUPLASRC)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(ALPINC) -I$(CUPLAINC) $(CUPLAFLAGS) -I$(BOOSTINC) -c $< -o $@

%.o : %.cc *.h ../../src/*.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(ALPINC) -I$(CUPLAINC) $(CUPLAFLAGS) -I$(BOOSTINC) -c $< -o $@

$(cxx_main): $(cxx_objects) check.o
	$(CXX) -o $@ $(cxx_objects) check.o -ldl $(LIBFLAGS)

.PHONY: clean

clean:
	cd ../../src && make clean
	rm -f *.o *.exe cupla/*.o cupla/manager/*.o

memcheck: $(cu_main)
	/usr/local/cuda/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes ./gcheck.exe 2 256 2

perf: force
	make clean && make
	time ./gcheck.exe -p 2048 256 12 && date

test: force
	./gcheck.exe -v 1 256 1

force:

#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
