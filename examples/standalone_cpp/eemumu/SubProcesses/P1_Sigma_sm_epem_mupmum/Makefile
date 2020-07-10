# standalone_cpp
LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_sm
CXXFLAGS= -g -O -I$(INCDIR) -I.
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

main=check.exe

all:
	cd ../../src && make
	make check

check: $(main)

objects=CPPProcess.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	cd ../../src && make clean
	rm -f $(main)
	rm -f $(objects)
