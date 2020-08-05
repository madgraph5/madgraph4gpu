#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// Memory layout for momenta
#define MGONGPU_LAYOUT_ASA 1 // default
//#define MGONGPU_LAYOUT_SOA 1
//#define MGONGPU_LAYOUT_AOS 1

// Curand random number generation
#define MGONGPU_CURAND_ONDEVICE 1 // default
//#define MGONGPU_CURAND_ONHOST 1

#ifdef __CUDACC__
#include <thrust/complex.h>
#else
#include <complex>
#endif

namespace mgOnGpu
{
  // Number of Events Per Page in the AOSOA (ASA) structure
  const int nepp = 32; // choose 32, like the number of threads in a warp

  // Maximum number of threads per block
  // ** NB kernels will statically allocate shared memory for all these threads
  const int ntpbMAX = 256;

  // Complex type
#ifdef __CUDACC__
  typedef thrust::complex<double> dcomplex;
#else
  typedef std::complex<double> dcomplex;
#endif

}

#endif // MGONGPUCONFIG_H
