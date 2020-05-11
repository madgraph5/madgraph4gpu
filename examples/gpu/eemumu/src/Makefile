CUARCHNUM=70
LIBDIR=../lib
CXXFLAGS= -O3 -I.
CUARCHFLAGS= -arch=compute_$(CUARCHNUM)
CUFLAGS= $(CUARCHFLAGS) -use_fast_math -rdc=true -lineinfo
CUINC=/usr/local/cuda/targets/x86_64-linux/include
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++

target=$(LIBDIR)/libmodel_sm.a
cxx_objects=Parameters_sm.o read_slha.o rambo.o
cu_objects=HelAmps_sm.o

all: $(target)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CUFLAGS:=$(filter-out -lineinfo,$(CUFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG
debug: CUFLAGS += -G
debug: $(target)

# sr fixme #
# not sure including the cuda includes here is a good idea...
# ... needed for thrust/complex which in principal is a C++ class
%.o : %.cc %.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(CUINC) -c $< -o $@

%.o : %.cu %.h
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) $(CUFLAGS) -c $< -o $@


$(target): $(cxx_objects) $(cu_objects)
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(cxx_objects) $(cu_objects)
