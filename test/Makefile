CXXFLAGS += -Igoogletest/googletest/include/ -std=c++11

all: googletest/build/lib/libgtest.a

googletest:
	git clone https://github.com/google/googletest.git -b release-1.11.0 googletest

googletest/build: googletest
	mkdir -p $@
	cd googletest/build && cmake -DBUILD_GMOCK=OFF ../

googletest/build/lib/libgtest.a: googletest/build
	$(MAKE) -C googletest/build

clean:
	rm -rf googletest
	rm -f src/*.o
