
CXX=g++

CXXFLAGS=-O3 -ffast-math --std=c++17 --openmp
LDFLAGS=--openmp

LIBDIR=../lib
target=$(LIBDIR)/libmodel_sa.a
cxx_objects=Parameters_sm.o read_slha.o

all: $(target)

%.o : %.cc %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(target): $(cxx_objects) $(cu_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target) *.o
