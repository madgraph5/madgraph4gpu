#ifndef MG5AMC_GPURUNTIME_H
#define MG5AMC_GPURUNTIME_H 1

// MG5AMC on GPU uses the CUDA runtime API, not the lower level CUDA driver API
// See https://docs.nvidia.com/cuda/cuda-runtime-api/driver-vs-runtime-api.html#driver-vs-runtime-api

#include "GpuAbstraction.h"

#include <iostream>

//--------------------------------------------------------------------------

// See https://stackoverflow.com/a/14038590
#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
#define checkGpu( code ) { assertGpu( code, __FILE__, __LINE__ ); }
inline void assertGpu( gpuError_t code, const char* file, int line, bool abort = true )
{
  if( code != gpuSuccess )
  {
    printf( "ERROR! assertGpu: '%s' (%d) in %s:%d\n", gpuGetErrorString( code ), code, file, line );
    if( abort ) assert( code == gpuSuccess );
  }
}
#endif /* clang-format on */

//--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
{
  // Instantiate a GpuRuntime at the beginnining of the application's main to
  // invoke gpuSetDevice(0) in the constructor and book a gpuDeviceReset() call in the destructor
  // *** FIXME! This will all need to be designed differently when going to multi-GPU nodes! ***
  struct GpuRuntime final
  {
    GpuRuntime( const bool debug = true )
      : m_debug( debug ) { setUp( m_debug ); }
    ~GpuRuntime() { tearDown( m_debug ); }
    GpuRuntime( const GpuRuntime& ) = delete;
    GpuRuntime( GpuRuntime&& ) = delete;
    GpuRuntime& operator=( const GpuRuntime& ) = delete;
    GpuRuntime& operator=( GpuRuntime&& ) = delete;
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
      if( debug ) std::cout << "__GpuRuntime: calling GpuSetDevice(0)" << std::endl;
      checkGpu( gpuSetDevice( 0 ) ); // SLOW!
    }

    // Tear down CUDA application (call cudaDeviceReset)
    // ** NB: strictly speaking this is not needed when using the CUDA runtime API **
    // Calling cudaDeviceReset on shutdown is only needed for checking memory leaks in cuda-memcheck
    // See https://docs.nvidia.com/cuda/cuda-memcheck/index.html#leak-checking
    static void tearDown( const bool debug = true )
    {
      if( debug ) std::cout << "__GpuRuntime: calling GpuDeviceReset()" << std::endl;
      checkGpu( gpuDeviceReset() );
    }
  };
}
#endif

//--------------------------------------------------------------------------

#endif // MG5AMC_GPURUNTIME_H
