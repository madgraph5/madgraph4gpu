//==========================================================================
// This file has been automatically generated for SYCL/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

//--------------------------------------------------------------------------

namespace Proc
{

%(process_class_definitions)s
  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          bool * isGoodHel,         // output: isGoodHel[ncomb] - device array
          size_t ievt,
          short *cHel,
          const fptype *cIPC,
          const fptype *cIPD
          );

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_setGoodHel(
          const bool * isGoodHel,  // input: isGoodHel[ncomb] - host array
          int * cNGoodHel_ptr,
          int* cGoodHel_ptr
          );

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          size_t ievt,
          short *cHel,
          const fptype *cIPC,
          const fptype *cIPD,
          int *cNGoodHel,
          int *cGoodHel
          );

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
