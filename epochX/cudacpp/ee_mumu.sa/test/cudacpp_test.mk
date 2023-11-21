# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Hageboeck (Dec 2020) for the CUDACPP plugin.
# Modified by: A. Valassi (2020-2023) for the CUDACPP plugin.

THISDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Compiler-specific googletest build directory (#125 and #738)
# In epochX, CXXNAMESUFFIX=_$(CXXNAME) is exported from cudacpp.mk
# In epoch1/epoch2, CXXNAMESUFFIX is undefined
$(info CXXNAMESUFFIX=$(CXXNAMESUFFIX))
BUILDDIR = build$(CXXNAMESUFFIX)
###$(info BUILDDIR=$(BUILDDIR))
INSTALLDIR = install$(CXXNAMESUFFIX)
###$(info INSTALLDIR=$(INSTALLDIR))

CXXFLAGS += -Igoogletest/googletest/include/ -std=c++11

all: googletest/$(INSTALLDIR)/lib64/libgtest.a

googletest/CMakeLists.txt:
	git clone https://github.com/google/googletest.git -b release-1.11.0 googletest

googletest/$(BUILDDIR)/Makefile: googletest/CMakeLists.txt
	mkdir -p googletest/$(BUILDDIR)
	cd googletest/$(BUILDDIR) && cmake -DCMAKE_INSTALL_PREFIX:PATH=$(THISDIR)/googletest/install -DBUILD_GMOCK=OFF ../

googletest/$(BUILDDIR)/lib/libgtest.a: googletest/$(BUILDDIR)/Makefile
	$(MAKE) -C googletest/$(BUILDDIR)

# NB 'make install' is no longer supported in googletest (issue 328)
# NB keep 'lib64' instead of 'lib' as in LCG cvmfs installations
googletest/$(INSTALLDIR)/lib64/libgtest.a: googletest/$(BUILDDIR)/lib/libgtest.a
	mkdir -p googletest/$(INSTALLDIR)/lib64
	cp googletest/$(BUILDDIR)/lib/lib*.a googletest/$(INSTALLDIR)/lib64/
	mkdir -p googletest/$(INSTALLDIR)/include
	cp -r googletest/googletest/include/gtest googletest/$(INSTALLDIR)/include/

clean:
	rm -rf googletest
