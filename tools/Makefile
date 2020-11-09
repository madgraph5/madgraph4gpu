# Set up tools, e.g. googletest.

.PHONY: gtest

all: gtest

googletest:
	git clone https://github.com/google/googletest.git -b release-1.10.0 googletest

googletest/build: googletest
	mkdir -p $@
	cd googletest/build && cmake -DBUILD_GMOCK=OFF ../

gtest: googletest/build
	$(MAKE) -C googletest/build

clean:
	rm -rf googletest
