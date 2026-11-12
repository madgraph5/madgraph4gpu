# Copyright (C) 2020-2025 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: S. Hageboeck (Dec 2020) for the CUDACPP plugin.
# Further modified by: S. Roiser, A. Valassi (2020-2025) for the CUDACPP plugin.
# Further modified by: A. Valassi (2026) for the MadVfold project (arXiv:2609.23500).

THISDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Host detection
UNAME_S := $(shell uname -s)
UNAME_P := $(shell uname -p)
UNAME_M := $(shell uname -m)

# Only add AVX2/FMA on non-mac and non-ARM hosts
ifeq ($(UNAME_S),Darwin)
  GTEST_CMAKE_FLAGS :=
else ifeq ($(UNAME_M),aarch64)
  GTEST_CMAKE_FLAGS :=
else
  GTEST_CMAKE_FLAGS := -DCMAKE_CXX_FLAGS="-mavx2 -mfma"
endif

# Compiler-specific googletest build directory (#125 and #738)
# In epochX, CXXNAMESUFFIX=_$(CXXNAME) is exported from cudacpp.mk
# In epoch1/epoch2, CXXNAMESUFFIX is undefined
$(info CXXNAMESUFFIX=$(CXXNAMESUFFIX))
BUILDDIR = build$(CXXNAMESUFFIX)
###$(info BUILDDIR=$(BUILDDIR))
INSTALLDIR = install$(CXXNAMESUFFIX)
###$(info INSTALLDIR=$(INSTALLDIR))

# [AV vfold] googletest is cloned, built and installed in GTEST/googletest
# GTEST is created by cudacpp.mk (target gtestlibs) as a symlink to a shared googletest directory or as a local directory
GOOGLETEST = GTEST/googletest

CXXFLAGS += -I$(GOOGLETEST)/googletest/include/ -std=c++11

all: $(GOOGLETEST)/$(INSTALLDIR)/lib64/libgtest.a

$(GOOGLETEST)/CMakeLists.txt:
	git clone https://github.com/google/googletest.git -b v1.17.0 $(GOOGLETEST)

$(GOOGLETEST)/$(BUILDDIR)/Makefile: $(GOOGLETEST)/CMakeLists.txt
	mkdir -p $(GOOGLETEST)/$(BUILDDIR)
	cd $(GOOGLETEST)/$(BUILDDIR) && cmake -DCMAKE_INSTALL_PREFIX:PATH=$(THISDIR)$(GOOGLETEST)/install $(GTEST_CMAKE_FLAGS) -DBUILD_GMOCK=OFF ../

$(GOOGLETEST)/$(BUILDDIR)/lib/libgtest.a: $(GOOGLETEST)/$(BUILDDIR)/Makefile
	$(MAKE) -C $(GOOGLETEST)/$(BUILDDIR)

# NB 'make install' is no longer supported in googletest (issue 328)
# NB keep 'lib64' instead of 'lib' as in LCG cvmfs installations
$(GOOGLETEST)/$(INSTALLDIR)/lib64/libgtest.a: $(GOOGLETEST)/$(BUILDDIR)/lib/libgtest.a
	mkdir -p $(GOOGLETEST)/$(INSTALLDIR)/lib64
	cp $(GOOGLETEST)/$(BUILDDIR)/lib/lib*.a $(GOOGLETEST)/$(INSTALLDIR)/lib64/
	mkdir -p $(GOOGLETEST)/$(INSTALLDIR)/include
	cp -r $(GOOGLETEST)/googletest/include/gtest $(GOOGLETEST)/$(INSTALLDIR)/include/

# [AV vfold] if GTEST is a symlink, only remove the symlink (keep the shared googletest installation)
clean:
	if [ -L GTEST ]; then rm -f GTEST; else rm -rf GTEST; fi
