//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#include "../../src/HelAmps_sm.cc"

#endif

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "CPPProcess.h"
#include "extras.h"

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

%(process_function_definitions)s
