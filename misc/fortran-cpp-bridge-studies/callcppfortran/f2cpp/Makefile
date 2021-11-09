LIBDIR=/opt/local/lib/gcc10

all:
	gfortran -c -g mg5_f.f
	g++ -c -g mg5_cpp.cpp
	g++ -o test mg5_f.o mg5_cpp.o -L${LIBDIR} -lgfortran

clean:
	rm -f *.o *.mod *.exe test
