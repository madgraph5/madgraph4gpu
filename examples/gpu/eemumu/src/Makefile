LIBDIR=../lib
CXXFLAGS= -g -O -I.
CUFLAGS= -arch=compute_35 -rdc=true
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++

target=$(LIBDIR)/libmodel_sm.a

all: $(target)

cxx_objects=Parameters_sm.o rambo.o read_slha.o
cu_objects=HelAmps_sm.o

%.o : %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o : %.cu
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@


$(target): $(cxx_objects) $(cu_objects)
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(cxx_objects) $(cu_objects)
