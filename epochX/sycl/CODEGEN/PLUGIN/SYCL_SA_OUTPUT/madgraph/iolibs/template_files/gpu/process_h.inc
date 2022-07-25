//==========================================================================
// This file has been automatically generated for SYCL standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

#include <vector>

//--------------------------------------------------------------------------

namespace Proc
{

%(process_class_definitions)s
  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,          // output: isGoodHel[ncomb] - device array
                            const short* __restrict__ cHel,
                            const fptype* __restrict__ cIPC,
                            const fptype* __restrict__ cIPD
                            );

  //--------------------------------------------------------------------------

  int sigmaKin_setGoodHel( const bool* isGoodHel, int* goodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  fptype sigmaKin( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                   const short* __restrict__ cHel,
                   const fptype* __restrict__ cIPC,
                   const fptype* __restrict__ cIPD,
                   const int* __restrict__ cNGoodHel,
                   const int* __restrict__ cGoodHel
                 );

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
