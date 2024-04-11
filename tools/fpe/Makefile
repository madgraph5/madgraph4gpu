all: testCC

%.o : %.cc
	$(CXX) -c $< -o $@

testCC : mainCC.o fpes.o
	$(CXX) -o $@ $^

clean:
	rm -f *.o testCC

