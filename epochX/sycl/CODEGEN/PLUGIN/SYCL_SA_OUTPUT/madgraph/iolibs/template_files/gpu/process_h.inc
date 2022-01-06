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

#ifdef SYCL_LANGUAGE_VERSION
  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool * isGoodHel,         // output: isGoodHel[ncomb] - device array
                            sycl::nd_item<3> item_ct1,
                            int *cHel,
                            const fptype *cIPC,
                            const fptype *cIPD
                            );
#else
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool* isGoodHel,           // output: isGoodHel[ncomb] - device array
                            const int nevt          // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
                            );
#endif

  //--------------------------------------------------------------------------

#ifdef SYCL_LANGUAGE_VERSION
  SYCL_EXTERNAL
  void sigmaKin_setGoodHel(const bool * isGoodHel, int * cNGoodHel_ptr, int* cGoodHel_ptr);  // input: isGoodHel[ncomb] - host array
#else
  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array
#endif

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs            // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                 #ifdef SYCL_LANGUAGE_VERSION
                     , sycl::nd_item<3> item_ct1,
                     int *cHel,
                     const fptype *cIPC,
                     const fptype *cIPD,
                     int *cNGoodHel,
                     int *cGoodHel
                 #else
                 , const int nevt  // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
                 #endif
                 );

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
