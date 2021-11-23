all:
	g++ -c mg5_cpp.cpp
	gfortran -c mg5_f90.f90
	gfortran -c mg5.f90
	gfortran -o mg5.exe mg5_cpp.o mg5_f90.o mg5.o -lstdc++

clean:
	rm -f *.o *.mod *.exe
