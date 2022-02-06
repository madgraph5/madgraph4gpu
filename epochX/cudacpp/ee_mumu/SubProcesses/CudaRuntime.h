#ifndef MG5AMC_CUDARUNTIME_H 
#define MG5AMC_CUDARUNTIME_H 1

// MG5AMC on GPU uses the CUDA runtime API, not the lower level CUDA driver API
// See https://docs.nvidia.com/cuda/cuda-runtime-api/driver-vs-runtime-api.html#driver-vs-runtime-api

#include <cassert>
#include <iostream>

//--------------------------------------------------------------------------

#ifdef __CUDACC__
  
#define checkCuda( code ) { assertCuda( code, __FILE__, __LINE__ ); }

inline void assertCuda( cudaError_t code, const char* file, int line, bool abort = true )
{
  if ( code != cudaSuccess )
  {
    printf( "ERROR! GPUassert: %s %s:%d\n", cudaGetErrorString(code), file, line );
    if ( abort ) assert( code == cudaSuccess );
  }
}

#endif

//--------------------------------------------------------------------------

#ifdef __CUDACC__
namespace mg5amcGpu
{

  // Book a cudaDeviceReset when CudaTearDown goes out of scope
  // This is needed to check for memory leaks in cuda-memcheck
  // See https://docs.nvidia.com/cuda/cuda-memcheck/index.html#leak-checking
  struct CudaTearDown {
    CudaTearDown(bool print) : _print(print) { }
    ~CudaTearDown() {
      //if ( _print ) std::cout << "Calling cudaDeviceReset()." << std::endl;
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
    bool _print{false};
  };

}
#endif

//--------------------------------------------------------------------------

#endif // MG5AMC_CUDARUNTIME_H
