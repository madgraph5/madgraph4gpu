LIBDIR=../lib
CXXFLAGS= -O -I.

target=$(LIBDIR)/libmodel_sm_ckm_no_b_mass.a

all: $(target)

objects=HelAmps_sm_ckm_no_b_mass.o Parameters_sm_ckm_no_b_mass.o rambo.o read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

