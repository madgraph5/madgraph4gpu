all: testCC testF

%.o : %.f
	$(FC) -c $< -o $@

%.o : %.cc
	$(CXX) -c $< -o $@

testCC : mainCC.o fpes.o
	$(CXX) -o $@ $^

testF : mainF.o fpes.o
	$(FC) -o $@ $^ -lstdc++

clean:
	rm -f *.o testCC testF

