all: testCC testF

# See https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html for the -ffpe options
# Default is -ffpe-summary=invalid,zero,overflow,underflow,denormal (all but inexact)
# Default is -ffpe-trap=none (none of invalid,zero,overflow,underflow,denormal,inexact)
%.o : %.f
	$(FC) -c $< -o $@ -ffpe-summary=invalid,zero,overflow #-ffpe-trap=invalid,zero,overflow

%.o : %.cc
	$(CXX) -c $< -o $@

testCC : mainCC.o fpes.o
	$(CXX) -o $@ $^

testF : mainF.o fpes.o
	$(FC) -o $@ $^ -lstdc++

clean:
	rm -f *.o testCC testF

