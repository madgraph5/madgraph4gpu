#ifndef MGKOKKOSTYPES_H
#define MGKOKKOSTYPES_H 1

#include "mgKokkosConfig.h"


#ifdef THRUST_COMPLEX
  #include <thrust/complex.h>
  template<typename T>
  using complex_t = thrust::complex<T>;
#else
  #include "Kokkos_Complex.hpp"
  template<typename T>
  using complex_t = Kokkos::complex<T>;
#endif

#ifndef __CUDACC__
#include <cmath>
using std::min;
using std::max;
using std::sqrt;
#endif

#ifdef __CUDACC__
#include <nvToolsExt.h> 
#else

inline void nvtxRangePush(const char* text){
  return;
}

inline void nvtxRangePop(void){
  return;
}
#endif // __CUDACC__


#endif // MGKOKKOSTYPES_H
