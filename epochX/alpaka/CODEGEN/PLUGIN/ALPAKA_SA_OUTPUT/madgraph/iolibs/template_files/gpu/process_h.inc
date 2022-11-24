//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
%(info_lines)s
//==========================================================================

#ifndef MG5_%(process_file_name)s_H
#define MG5_%(process_file_name)s_H

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#ifdef ALPAKA
#include <cuda_to_cupla.hpp>
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "Parameters_%(model_name)s.h"

//--------------------------------------------------------------------------

#if defined(ALPAKA) || defined(__CUDACC__)

#define checkCuda( code ) { assertCuda( code, __FILE__, __LINE__ ); }

inline void assertCuda( cudaError_t code, const char* file, int line, bool abort = true )
{
  if ( code != cudaSuccess )
  {
    printf( "GPUassert: %%s %%s:%%d\n", cudaGetErrorString(code), file, line );
    if ( abort ) assert( code == cudaSuccess );
  }
}

#endif

//--------------------------------------------------------------------------

#ifdef ALPAKA
namespace gProc
{
%(process_class_definitions)s
  //--------------------------------------------------------------------------

  struct sigmaKin_getGoodHel
  {
    template< typename T_Acc >
    ALPAKA_FN_ACC
    void operator()( T_Acc const &acc,
                     const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                     fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                     bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                   ) const;
  };

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  struct sigmaKin
  {
    template< typename T_Acc >
    ALPAKA_FN_ACC
    void operator()( T_Acc const &acc,
                     const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                     fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                   ) const;
  };

  //--------------------------------------------------------------------------
}
#else
#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{

%(process_class_definitions)s
  //--------------------------------------------------------------------------

  __global__
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
#ifndef __CUDACC__
                            , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                            );

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ); // input: isGoodHel[ncomb] - host array

  //--------------------------------------------------------------------------

  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 );

  //--------------------------------------------------------------------------
}
#endif

#endif // MG5_%(process_file_name)s_H
