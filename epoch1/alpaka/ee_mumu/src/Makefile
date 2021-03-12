LIBDIR=../lib
CXXFLAGS= -O3 -std=c++14 -I. -Wall -Wshadow 
CXX=g++

target=$(LIBDIR)/libmodel_sm.a
cxx_objects=Parameters_sm.o read_slha.o alsimple.o

all: $(target)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -g -O0 -DDEBUG2
debug: $(target)

# NB: cuda includes are needed in the C++ code for curand.h
%.o : %.cc *.h
	$(CXX) -I$(ALPINC) -I$(BOOSTINC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(target): $(cxx_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(cxx_objects)
