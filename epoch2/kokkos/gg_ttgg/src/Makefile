
CXX=gcc
ICPX=icpx

CXXFLAGS=-O3 -ffast-math --std=c++17

LIBDIR=../lib
gen_lib=$(LIBDIR)/libmodel_sa.a
iris_lib=$(LIBDIR)/libmodel_sa_iris.a
cxx_objects=Parameters_sm.o read_slha.o
iris_objs=iParameters_sm.o iread_slha.o

all: generic

generic: $(gen_lib)

iris: $(iris_lib)

%.o : %.cc %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(gen_lib): $(cxx_objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(cxx_objects) $(cu_objects)
	ranlib $@

iParameters_sm.o: Parameters_sm.cc Parameters_sm.h
	$(ICPX) $(CXXFLAGS) -c $< -o $@

iread_slha.o: read_slha.cc read_slha.h
	$(ICPX) $(CXXFLAGS) -c $< -o $@


$(iris_lib): iParameters_sm.o iread_slha.o
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

.PHONY: clean

clean:
	rm -f $(target) *.o
