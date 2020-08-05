CUARCHNUM=70
LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -O3 -I$(INCDIR) -I.
CUARCHFLAGS= -arch=compute_$(CUARCHNUM) #--gpu-architecture=compute_35 --gpu-code=sm_35
CUFLAGS= $(CUARCHFLAGS) -use_fast_math -lineinfo -DUSE_NVTX --default-stream per-thread
CUINC=/usr/local/cuda/targets/x86_64-linux/include
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)
CULIBFLAGS= -L/usr/local/cuda/targets/x86_64-linux/lib
NVCC=nvcc
CXX=g++
MAKEDEBUG=

cu_main=gcheck.exe
cxx_main=check.exe

cu_objects=gCPPProcess.o gcheck.o
cxx_objects=CPPProcess.o

# Assuming uname is available, detect if architecture is power
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
    CUFLAGS+= -Xcompiler -mno-float128
endif

all:
	cd ../../src && make
	make $(cu_main)
	make $(cxx_main)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG2
debug: CUFLAGS:=$(filter-out -lineinfo,$(CUFLAGS))
debug: CUFLAGS += -G
debug: MAKEDEBUG := debug
debug: $(cu_main)

$(LIBDIR)/lib$(MODELLIB).a:
	@cd ../../src && make $(MAKEDEBUG)

gcheck.o: gcheck.cu *.h ../../src/*.h ../../src/*.cu
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@

%.o : %.cu *.h ../../src/*.h
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@

%.o : %.cc *.h ../../src/*.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CUINC) -c $< -o $@

$(cu_main): $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) $(cu_objects)
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) $(LIBFLAGS) $(CUARCHFLAGS) -lcuda -lcurand -lcudart

$(cxx_main): $(LIBDIR)/lib$(MODELLIB).a $(cxx_objects) check.o
	$(CXX) -o $@ $(cxx_objects) check.o $(LIBFLAGS) $(CULIBFLAGS) -lcurand

.PHONY: clean

clean:
	cd ../../src && make clean
	rm -f $(cu_main)
	rm -f $(cu_objects) $(cxx_objects)

memcheck: $(cu_main)
	/usr/local/cuda/bin/cuda-memcheck --check-api-memory-access yes --check-deprecated-instr yes --check-device-heap yes --demangle full --language c --leak-check full --racecheck-report all --report-api-errors all --show-backtrace yes --tool memcheck --track-unused-memory yes ./gcheck.exe 2 32 2

perf: force
	make clean && make
	time ./gcheck.exe -p 2048 256 12 && date

test: force
	./gcheck.exe -v 1 32 1

force:

#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
