#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// Memory layout for momenta
//#define MGONGPU_LAYOUT_ASA 1 // default
//#define MGONGPU_LAYOUT_SOA 1
#define MGONGPU_LAYOUT_AOS 1

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
  // --- Physics process-specific constants that are best declared at compile time

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)

  const int npari = 2; // #particles in the initial state (incoming): e+ e-
  const int nparf = 2; // #particles in the final state (outgoing): mu+ mu-
  const int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  // --- Platform-specific software implementation details

  // Maximum number of threads per block
  // ** NB kernels will statically allocate shared memory for all these threads
  const int ntpbMAX = 256;

  // Number of Events Per Page in the AOSOA (ASA) structure
  //const int nepp = 32; // choose 32, i.e. the number of threads in a warp
  const int nepp = ntpbMAX; // choose 256, i.e. the max number of threads in a block

  // Complex type
#ifdef __CUDACC__
  typedef thrust::complex<double> dcomplex; // two doubles: RI
#else
  typedef std::complex<double> dcomplex; // two doubles: RI
#endif

  // Vector types: <type>_v is a <type>[256]
  typedef double double_v[ntpbMAX];
  typedef dcomplex dcomplex_v[ntpbMAX]; // RIRIRIRI: eventually move to RRRRIIII?

}

#endif // MGONGPUCONFIG_H
