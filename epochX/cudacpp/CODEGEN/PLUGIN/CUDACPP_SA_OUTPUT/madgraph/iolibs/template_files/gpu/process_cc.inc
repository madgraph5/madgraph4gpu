//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"

#include "CPPProcess.h"
#include "CudaRuntime.h"
%(hel_amps_h)s
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

%(process_function_definitions)s
