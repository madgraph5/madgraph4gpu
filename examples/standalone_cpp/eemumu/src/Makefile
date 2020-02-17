LIBDIR=../lib
CXXFLAGS= -O -I.

target=$(LIBDIR)/libmodel_sm.a

all: $(target)

objects=HelAmps_sm.o Parameters_sm.o rambo.o read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

