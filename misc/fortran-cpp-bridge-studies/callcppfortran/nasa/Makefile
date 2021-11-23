all:
	g++ -c pointrectangle_class.cc
	gfortran -c pointrectangle_module.f90
	gfortran -c pointrectangle_main.f90
	gfortran -o pointrectangle.ex pointrectangle_class.o pointrectangle_module.o pointrectangle_main.o -lstdc++

clean:
	rm -f *.o *.ex *.mod
