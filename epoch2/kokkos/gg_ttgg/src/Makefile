
CXX=gcc
ICPX=icpx
CLANG=clang++
HIPCC=hipcc

CXXFLAGS=-O3 -ffast-math --std=c++17 -lineinfo

LIBDIR=../lib
gen_lib=$(LIBDIR)/libmodel_sa.a
intel_lib=$(LIBDIR)/libmodel_sa_intel.a
hip_lib=$(LIBDIR)/libmodel_sa_hip.a
cxx_objects=Parameters_sm.o read_slha.o
intel_objs=iParameters_sm.o iread_slha.o
hip_objs=hParameters_sm.o hread_slha.o

all: generic

generic: $(gen_lib)

intel: $(intel_lib)

hip: $(hip_lib)

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


$(intel_lib): iParameters_sm.o iread_slha.o
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@

hParameters_sm.o: Parameters_sm.cc Parameters_sm.h
	$(HIPCC) $(CXXFLAGS) -c $< -o $@

hread_slha.o: read_slha.cc read_slha.h
	$(HIPCC) $(CXXFLAGS) -c $< -o $@


$(hip_lib): hParameters_sm.o hread_slha.o
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $^
	ranlib $@



.PHONY: clean

clean:
	rm -f $(target) *.o
