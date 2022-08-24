//==========================================================================
// This file has been automatically generated for SYCL standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

#include <vector>

//--------------------------------------------------------------------------

namespace mg5amcGpu
{

%(process_class_definitions)s
  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta[nevt*npar*4]
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool* isGoodHel,          // output: isGoodHel[ncomb] - device array
                            const size_t ievt,
                            short* cHel,
                            const fptype* cIPC,
                            const fptype* cIPD
                            );

  //--------------------------------------------------------------------------

  int sigmaKin_setGoodHel( const bool* isGoodHel, int* goodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin( const fptype* allmomenta, // input: momenta[nevt*npar*4]
                 fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                 size_t ievt,
                 short* cHel,
                 const fptype* cIPC,
                 const fptype* cIPD,
                 int* cNGoodHel,
                 int* cGoodHel
                 );

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
