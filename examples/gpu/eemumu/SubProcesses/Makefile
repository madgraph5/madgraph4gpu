CUARCHNUM=70
LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -O3 -I$(INCDIR) -I.
CUARCHFLAGS= -arch=compute_$(CUARCHNUM) #--gpu-architecture=compute_35 --gpu-code=sm_35
CUFLAGS= $(CUARCHFLAGS) -rdc=true
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++
MAKEDEBUG=

main=check.exe

all:
	cd ../../src && make
	make check

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -g -O0
debug: CUFLAGS += -G
debug: MAKEDEBUG := debug
debug: check

check: $(main)

cu_objects=CPPProcess.o check_sa.o

%.o : %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o : %.cu
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@

$(main): $(cxx_objects) $(cu_objects) $(LIBDIR)/lib$(MODELLIB).a
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) $(LIBFLAGS) $(CUARCHFLAGS) -lcuda

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src && make $(MAKEDEBUG)

.PHONY: clean

clean:
	cd ../../src && make clean
	rm -f $(main)
	rm -f $(cu_objects) $(cxx_objects)


#Allowed values for this option: 'compute_30', 'compute_32', 'compute_35', 'compute_37', 'compute_50', 'compute_52', 'compute_53', 'compute_60', 'compute_61', 'compute_62', 'compute_70', 'compute_72', 'compute_75', 'sm_30', 'sm_32', 'sm_35', 'sm_37', 'sm_50', 'sm_52', 'sm_53', 'sm_60', 'sm_61', 'sm_62', 'sm_70', 'sm_72', 'sm_75'.

# Max compute architectures
# cern batch (tesla v100): 70
# jetson nano (maxwell): 35
