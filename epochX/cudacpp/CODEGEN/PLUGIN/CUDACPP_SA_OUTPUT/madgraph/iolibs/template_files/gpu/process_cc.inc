//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "CudaRuntime.h"
%(hel_amps_h)s
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"
#include "Parameters_sm.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

%(process_function_definitions)s
//==========================================================================
