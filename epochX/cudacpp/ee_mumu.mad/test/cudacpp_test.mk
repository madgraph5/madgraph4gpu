# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Hageboeck (Dec 2020) for the CUDACPP plugin.
# Modified by: A. Valassi (2020-2023) for the CUDACPP plugin.

THISDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CXXFLAGS += -Igoogletest/googletest/include/ -std=c++11

all: googletest/build/lib/libgtest.a

googletest:
	git clone https://github.com/google/googletest.git -b release-1.11.0 googletest

googletest/build: googletest
	mkdir -p $@
	cd googletest/build && cmake -DCMAKE_INSTALL_PREFIX:PATH=$(THISDIR)/googletest/install -DBUILD_GMOCK=OFF ../

googletest/build/lib/libgtest.a: googletest/build
	$(MAKE) -C googletest/build
	$(MAKE) -C googletest/build install

clean:
	rm -rf googletest
