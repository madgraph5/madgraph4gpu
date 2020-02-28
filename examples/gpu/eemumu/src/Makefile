LIBDIR=../lib
CXXFLAGS= -g -O0 -I.
CUFLAGS= -arch=compute_35 -rdc=true -G
# CUFLAGS= -arch=compute_35 -G
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++

target=$(LIBDIR)/libmodel_sm.a

all: $(target)

cxx_objects=Parameters_sm.o read_slha.o rambo.o
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
