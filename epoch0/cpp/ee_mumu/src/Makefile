LIBDIR=../lib
CXXFLAGS= -g -O -I.

target=$(LIBDIR)/libmodel_sm.a

all: $(target)

objects=HelAmps_sm.o Parameters_sm.o rambo.o read_slha.o

$(target): $(objects)
	if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR); fi
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

