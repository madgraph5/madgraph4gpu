CUARCHNUM=75
LIBDIR=../lib
CXXFLAGS= -O3 -I. -I/home/jchilders/git/kokkos/install/include --openmp --std=c++17
CUARCHFLAGS= -arch=compute_$(CUARCHNUM) --expt-extended-lambda
CUFLAGS= $(CUARCHFLAGS) -use_fast_math -lineinfo
CUINC=/usr/local/cuda/targets/x86_64-linux/include
NVCC=/home/jchilders/git/kokkos/bin/nvcc_wrapper
CXX=g++

# Assuming uname is available, detect if architecture is power
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),ppc64le)
    CUFLAGS+= -Xcompiler -mno-float128
endif


target=$(LIBDIR)/libmodel_sm.a
cxx_objects=Parameters_sm.o read_slha.o rambo.o
cu_objects=

all: $(target)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS)) -I/home/jchilders/git/kokkos/install/include --openmp --std=c++17
debug: CUFLAGS:=$(filter-out -lineinfo,$(CUFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG2
debug: CUFLAGS += -G
debug: $(target)

# sr fixme #
# not sure including the cuda includes here is a good idea...
# ... needed for thrust/complex which in principal is a C++ class
%.o : %.cc %.h
	$(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) -I$(CUINC)

%.o : %.cu %.h
	$(NVCC) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS)


$(target): $(cxx_objects) $(cu_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(cxx_objects) $(cu_objects)
