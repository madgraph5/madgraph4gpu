CXXFLAGS += -Igoogletest/googletest/include/ -std=c++11

all: gtest src/MadgraphTest.o

%.o: %.cc include/*.h
	$(CXX) $< -c -o $@ $(CPPFLAGS) $(CXXFLAGS) -Iinclude/

.PHONY: gtest

googletest:
	git clone https://github.com/google/googletest.git -b release-1.10.0 googletest

googletest/build: googletest
	mkdir -p $@
	cd googletest/build && cmake -DBUILD_GMOCK=OFF ../

gtest: googletest/build
	$(MAKE) -C googletest/build

clean:
	rm -rf googletest
	rm src/*.o
