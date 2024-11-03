# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Mar 2024) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#-------------------------------------------------------------------------------

#=== Check that the user-defined choices of BACKEND, FPTYPE, HELINL, HRDCOD are supported
#=== Configure default values for these variables if no user-defined choices exist

# Set the default BACKEND (CUDA, HIP or C++/SIMD) choice
ifeq ($(BACKEND),)
  override BACKEND = cppauto
endif

# Set the default FPTYPE (floating point type) choice
# NB: this only affects manual 'make' builds (madevent 'launch' builds are controlled by floating_type in run_card.dat)
ifeq ($(FPTYPE),)
  # OLD DEFAULT UP TO v1.00.00 INCLUDED (inconsistent with default floating_type='m' in run_card.dat)
  ###override FPTYPE = d
  # NEW DEFAULT (#995) AS OF v1.00.01 (now consistent with default floating_type='m' in run_card.dat)
  override FPTYPE = m
endif

# Set the default HELINL (inline helicities?) choice
ifeq ($(HELINL),)
  override HELINL = 0
endif

# Set the default HRDCOD (hardcode cIPD physics parameters?) choice
ifeq ($(HRDCOD),)
  override HRDCOD = 0
endif

# Check that the user-defined choices of BACKEND, FPTYPE, HELINL, HRDCOD are supported
# (NB: use 'filter' and 'words' instead of 'findstring' because they properly handle whitespace-separated words)
override SUPPORTED_BACKENDS = cuda hip cppnone cppsse4 cppavx2 cpp512y cpp512z cppauto
ifneq ($(words $(filter $(BACKEND), $(SUPPORTED_BACKENDS))),1)
  $(error Invalid backend BACKEND='$(BACKEND)': supported backends are $(foreach backend,$(SUPPORTED_BACKENDS),'$(backend)'))
endif

override SUPPORTED_FPTYPES = d f m
ifneq ($(words $(filter $(FPTYPE), $(SUPPORTED_FPTYPES))),1)
  $(error Invalid fptype FPTYPE='$(FPTYPE)': supported fptypes are $(foreach fptype,$(SUPPORTED_FPTYPES),'$(fptype)'))
endif

override SUPPORTED_HELINLS = 0 1
ifneq ($(words $(filter $(HELINL), $(SUPPORTED_HELINLS))),1)
  $(error Invalid helinl HELINL='$(HELINL)': supported helinls are $(foreach helinl,$(SUPPORTED_HELINLS),'$(helinl)'))
endif

override SUPPORTED_HRDCODS = 0 1
ifneq ($(words $(filter $(HRDCOD), $(SUPPORTED_HRDCODS))),1)
  $(error Invalid hrdcod HRDCOD='$(HRDCOD)': supported hrdcods are $(foreach hrdcod,$(SUPPORTED_HRDCODS),'$(hrdcod)'))
endif

# Print out BACKEND, FPTYPE, HELINL, HRDCOD
###$(info BACKEND='$(BACKEND)')
###$(info FPTYPE='$(FPTYPE)')
###$(info HELINL='$(HELINL)')
###$(info HRDCOD='$(HRDCOD)')

#-------------------------------------------------------------------------------

# Stop immediately if BACKEND=cuda but nvcc is missing
ifeq ($(BACKEND),cuda)
  ifeq ($(shell which nvcc 2>/dev/null),)
    $(error BACKEND=$(BACKEND) but nvcc was not found)
  endif
endif

# Stop immediately if BACKEND=hip but hipcc is missing
ifeq ($(BACKEND),hip)
  ifeq ($(shell which hipcc 2>/dev/null),)
    $(error BACKEND=$(BACKEND) but hipcc was not found)
  endif
endif

#-------------------------------------------------------------------------------

#=== Configure CUDACPP_BUILDDIR

# Build directory "short" tag (defines target and path to the optional build directory)
# (Rationale: keep directory names shorter, e.g. do not include random number generator choice)
# ** NB: using ':=' here ensures that 'cppauto' is used as such before being changed later on!
override DIRTAG := $(patsubst cpp%,%,$(BACKEND))_$(FPTYPE)_inl$(HELINL)_hrd$(HRDCOD)

# Build directory: current directory by default, or build.$(DIRTAG) if USEBUILDDIR==1
ifeq ($(USEBUILDDIR),1)
  override CUDACPP_BUILDDIR = build.$(DIRTAG)
else
  override CUDACPP_BUILDDIR = .
endif
###$(info USEBUILDDIR='$(USEBUILDDIR)')
###$(info CUDACPP_BUILDDIR='$(CUDACPP_BUILDDIR)')

#-------------------------------------------------------------------------------
