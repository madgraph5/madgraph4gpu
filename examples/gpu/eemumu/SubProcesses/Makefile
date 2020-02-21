# gpu
LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -g -O -I$(INCDIR) -I.
CUFLAGS= -arch=compute_35 -rdc=true
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++

main=check.exe

all:
	cd ../../src && make
	make check

check: $(main)

# cxx_objects=check_sa.o
cu_objects=CPPProcess.o check_sa.o

%.o : %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o : %.cu
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@

$(main): $(cxx_objects) $(cu_objects) $(LIBDIR)/lib$(MODELLIB).a
	$(NVCC) -o $@ $(cxx_objects) $(cu_objects) $(LIBFLAGS) -arch=compute_35

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	cd ../../src && make clean
	rm -f $(main)
	rm -f $(cu_objects) $(cxx_objects)
