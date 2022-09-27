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
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#endif

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

%(process_function_definitions)s


   __global__ normalise_output(fptype* allMEs, const fptype* allNumetators, cont fptype* allDenominators){

	    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
            const int globadenom = %(den_factors)s;
            //todo get ievt;
            allME[ievt] = allME[ievt] * allNumerators[ievt]/(allDenominators[ievt]*globaldenom);
            return
}

//==========================================================================
