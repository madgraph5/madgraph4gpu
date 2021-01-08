#CUARCHNUM=70
LIBDIR=../lib
#CXXFLAGS= -O3 -I.
CXXFLAGS=-Wno-sycl-strict -I.
#CUARCHFLAGS= -arch=compute_$(CUARCHNUM)
#CUFLAGS= $(CUARCHFLAGS) -use_fast_math -lineinfo
#CUINC=/usr/local/cuda/targets/x86_64-linux/include
CXX=dpcpp


target=$(LIBDIR)/libmodel_sm.a
cxx_objects=Parameters_sm.o read_slha.o rambo.o


all: $(target)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG2
debug: $(target)

# sr fixme #
# not sure including the cuda includes here is a good idea...
# ... needed for thrust/complex which in principal is a C++ class
%.o : %.cc %.h
	$(CXX)  $(CPPFLAGS) $(CXXFLAGS)  -c $< -o $@



$(target): $(cxx_objects) 
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) 
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(cxx_objects) 
