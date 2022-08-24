//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

#include <vector>

//--------------------------------------------------------------------------

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
%(process_class_definitions)s
  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__ void
  computeDependentCouplings( const fptype* allgs,    // input: Gs[nevt]
                             fptype* allcouplings ); // output: couplings[nevt*ndcoup*2]
#else
  __global__ void
  computeDependentCouplings( const fptype* allgs,  // input: Gs[nevt]
                             fptype* allcouplings, // output: couplings[nevt*ndcoup*2]
                             const int nevt );     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__ /* clang-format off */
  __global__ void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel );          // output: isGoodHel[ncomb] - device array
#else
  __global__ void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - device array
                       const int nevt );           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif /* clang-format on */

  //--------------------------------------------------------------------------

  void
  sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

#ifdef __CUDACC__ /* clang-format off */
  __global__ void
  sigmaKin( const fptype* allmomenta,       // input: momenta[nevt*npar*4]
            const fptype* allcouplings,     // input: couplings[nevt*ndcoup*2]
            fptype* allMEs                  // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            , fptype* allNumerators         // output: multichannel numerators[nevt], running_sum_over_helicities
            , fptype* allDenominators       // output: multichannel denominators[nevt], running_sum_over_helicities
            , const unsigned int channelId  // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
            );
#else
  __global__ void
  sigmaKin( const fptype* allmomenta,     // input: momenta[nevt*npar*4]
            const fptype* allcouplings,   // input: couplings[nevt*ndcoup*2]
            fptype* allMEs,               // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            fptype* allNumerators,        // output: multichannel numerators[nevt], running_sum_over_helicities
            fptype* allDenominators,      // output: multichannel denominators[nevt], running_sum_over_helicities
            const unsigned int channelId, // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
            const int nevt );             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif /* clang-format on */

  //--------------------------------------------------------------------------
}

#endif // MG5_%(process_file_name)s_H
