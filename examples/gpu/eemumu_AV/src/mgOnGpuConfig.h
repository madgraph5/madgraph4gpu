#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// Memory layout for momenta (CHOOSE ONLY ONE)
#define MGONGPU_LAYOUT_ASA 1 // default
//#define MGONGPU_LAYOUT_SOA 1
//#define MGONGPU_LAYOUT_AOS 1

// Curand random number generation (CHOOSE ONLY ONE)
#define MGONGPU_CURAND_ONDEVICE 1 // default
//#define MGONGPU_CURAND_ONHOST 1

// Use global memory or shared memory for wavefunctions (CHOOSE ONLY ONE)
#define MGONGPU_WFMEM_LOCAL 1 // default
//#define MGONGPU_WFMEM_GLOBAL 1
//#define MGONGPU_WFMEM_SHARED 1

// Floating point precision (CHOOSE ONLY ONE)
//#define MGONGPU_FPTYPE_DOUBLE 1 // default
#define MGONGPU_FPTYPE_FLOAT 1

#ifdef __CUDACC__
#include <thrust/complex.h>
#else
#include <complex>
#endif

namespace mgOnGpu
{
  // --- Physics process-specific constants that are best declared at compile time

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
  const int nw6 = 6; // dimension of each wavefunction (see KEK 91-11)

  const int npari = 2; // #particles in the initial state (incoming): e+ e-
  const int nparf = 2; // #particles in the final state (outgoing): mu+ mu-
  const int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  const int nwf = 5; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  const int nbpgMAX = 2048;

  // Maximum number of threads per block
  const int ntpbMAX = 256;
  //const int ntpbMAX = 8; // FOR DEBUGGING!

  // Number of Events Per Page in the AOSOA (ASA) structure
  // ** TODO: this will disappear and become dynamic...
  const int nepp = 32; // choose 32, i.e. the number of threads in a warp
  //const int nepp = ntpbMAX; // choose 256, i.e. the max number of threads in a block
  //const int nepp = 4; // FOR DEBUGGING!

  // Floating point type: fptype (tflpoint? tfloatpt?)
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef double fptype; // double precision (8 bytes, fp64)
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef float fptype; // single precision (4 bytes, fp32)
#endif

  // Complex type: cxtype (tcomplex?)
#ifdef __CUDACC__
  typedef thrust::complex<fptype> cxtype; // two doubles: RI
#else
  typedef std::complex<fptype> cxtype; // two doubles: RI
#endif

  // Vector types: <type>_v is a <type>[256]
  //typedef double double_v[ntpbMAX];
  //typedef cxtype cxtype_v[ntpbMAX]; // RIRIRIRI: eventually move to RRRRIIII?

}

// Expose typedefs outside the namespace
using mgOnGpu::fptype;
using mgOnGpu::cxtype;

#endif // MGONGPUCONFIG_H
