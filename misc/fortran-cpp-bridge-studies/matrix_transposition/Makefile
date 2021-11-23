# CXXFLAGS =

all: main

debug: CXXFLAGS += -DDEBUG -g
debug: main

main: matrix.cu
	nvcc $(CXXFLAGS) -o main $ matrix.cu

clean:
	rm -f main
