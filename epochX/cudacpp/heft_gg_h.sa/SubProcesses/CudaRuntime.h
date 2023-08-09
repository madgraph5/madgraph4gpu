// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Roiser (Jul 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.

#ifndef MG5AMC_CUDARUNTIME_H
#define MG5AMC_CUDARUNTIME_H 1

// MG5AMC on GPU uses the CUDA runtime API, not the lower level CUDA driver API
// See https://docs.nvidia.com/cuda/cuda-runtime-api/driver-vs-runtime-api.html#driver-vs-runtime-api

#include <cassert>
#include <iostream>

//--------------------------------------------------------------------------

// See https://stackoverflow.com/a/14038590
#ifdef __CUDACC__ /* clang-format off */
#define checkCuda( code ) { assertCuda( code, __FILE__, __LINE__ ); }
inline void assertCuda( cudaError_t code, const char* file, int line, bool abort = true )
{
  if( code != cudaSuccess )
  {
    printf( "ERROR! assertCuda: '%s' (%d) in %s:%d\n", cudaGetErrorString( code ), code, file, line );
    if( abort ) assert( code == cudaSuccess );
  }
}
#endif /* clang-format on */

//--------------------------------------------------------------------------

#ifdef __CUDACC__
namespace mg5amcGpu
{
  // Instantiate a CudaRuntime at the beginnining of the application's main to
  // invoke cudaSetDevice(0) in the constructor and book a cudaDeviceReset() call in the destructor
  // *** FIXME! This will all need to be designed differently when going to multi-GPU nodes! ***
  struct CudaRuntime final
  {
    CudaRuntime( const bool debug = true )
      : m_debug( debug ) { setUp( m_debug ); }
    ~CudaRuntime() { tearDown( m_debug ); }
    CudaRuntime( const CudaRuntime& ) = delete;
    CudaRuntime( CudaRuntime&& ) = delete;
    CudaRuntime& operator=( const CudaRuntime& ) = delete;
    CudaRuntime& operator=( CudaRuntime&& ) = delete;
    bool m_debug;

    // Set up CUDA application
    // ** NB: strictly speaking this is not needed when using the CUDA runtime API **
    // Calling cudaSetDevice on startup is useful to properly book-keep the time spent in CUDA initialization
    static void setUp( const bool debug = true )
    {
      // ** NB: it is useful to call cudaSetDevice, or cudaFree, to properly book-keep the time spent in CUDA initialization
      // ** NB: otherwise, the first CUDA operation (eg a cudaMemcpyToSymbol in CPPProcess ctor) appears to take much longer!
      /*
      // [We initially added cudaFree(0) to "ease profile analysis" only because it shows up as a big recognizable block!]
      // No explicit initialization is needed: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#initialization
      // It is not clear what cudaFree(0) does at all: https://stackoverflow.com/questions/69967813/
      if ( debug ) std::cout << "__CudaRuntime: calling cudaFree(0)" << std::endl;
      checkCuda( cudaFree( 0 ) ); // SLOW!
      */
      // Replace cudaFree(0) by cudaSetDevice(0), even if it is not really needed either
      // (but see https://developer.nvidia.com/blog/cuda-pro-tip-always-set-current-device-avoid-multithreading-bugs)
      if( debug ) std::cout << "__CudaRuntime: calling cudaSetDevice(0)" << std::endl;
      checkCuda( cudaSetDevice( 0 ) ); // SLOW!
    }

    // Tear down CUDA application (call cudaDeviceReset)
    // ** NB: strictly speaking this is not needed when using the CUDA runtime API **
    // Calling cudaDeviceReset on shutdown is only needed for checking memory leaks in cuda-memcheck
    // See https://docs.nvidia.com/cuda/cuda-memcheck/index.html#leak-checking
    static void tearDown( const bool debug = true )
    {
      if( debug ) std::cout << "__CudaRuntime: calling cudaDeviceReset()" << std::endl;
      checkCuda( cudaDeviceReset() );
    }
  };

}
#endif

//--------------------------------------------------------------------------

#endif // MG5AMC_CUDARUNTIME_H
