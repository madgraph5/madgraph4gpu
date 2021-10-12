#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// ** NB1 Throughputs (e.g. 6.8E8) are events/sec for "./gcheck.exe -p 65536 128 12"
// ** NB2 Baseline on b7g47n0004 fluctuates (probably depends on load on other VMs)

// Choose how random numbers are generated
// If one of these macros has been set from outside with e.g. -DMGONGPU_COMMONRAND_ONHOST, nothing happens
#if not defined MGONGPU_CURAND_ONDEVICE and not defined MGONGPU_CURAND_ONHOST and not defined MGONGPU_COMMONRAND_ONHOST
// Curand random number generation (CHOOSE ONLY ONE)
#define MGONGPU_CURAND_ONDEVICE 1 // default (curand: CUDA on device, C++ on host)
//#define MGONGPU_CURAND_ONHOST 1 // (curand: CUDA on host, C++ on host)
//#define MGONGPU_COMMONRAND_ONHOST 1 // (common rand: CUDA on host, C++ on host)
#endif

// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default (~6.8E8)
//#define MGONGPU_FPTYPE_FLOAT 1 // 2.4x faster (~1.64E9 against 6.8E8)

// Complex type in cuda: thrust or cucomplex (CHOOSE ONLY ONE)
#ifdef __CUDACC__
#define MGONGPU_CXTYPE_THRUST 1 // default (~6.8E8)
//#define MGONGPU_CXTYPE_CUCOMPLEX 1 // ~5 percent slower (6.5E8 against 6.8E8)
#endif

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
#ifdef __CUDACC__
#undef MGONGPU_NSIGHT_DEBUG // default
//#define MGONGPU_NSIGHT_DEBUG 1
#endif

namespace mgOnGpu
{

  // --- Type definitions

  // Floating point type: fptype
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef double fptype; // double precision (8 bytes, fp64)
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef float fptype; // single precision (4 bytes, fp32)
#endif

  // --- Physics process-specific constants that are best declared at compile time

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
  const int nw6 = %(wavefuncsize)d; // dimension of each wavefunction (see KEK 91-11)

  const int npari = %(nincoming)d; // #particles in the initial state (incoming): e+ e-
  const int nparf = %(noutcoming)d; // #particles in the final state (outgoing): mu+ mu-
  const int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  const int nwf = %(nwavefunc)d; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)

  const int ncomb = %(nbhel)d; // #helicity combinations: 16=2(spin up/down for fermions)**4(npar)

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  const int ntpbMAX = 256;

  // Vector sizes for AOSOA memory layouts (GPU coalesced memory access, CPU SIMD vectorization)
  // (these are all best kept as a compile-time constants: see issue #23)
  // There is no vectorization in epochX yet...
#undef MGONGPU_CPPSIMD

  // Number of Events Per Page in the momenta AOSOA (ASA) structure
  // (this is best kept as a compile-time constant: see issue #23)
#if defined MGONGPU_FPTYPE_DOUBLE
  const int neppM = 4; // DEFAULT: one 32-byte cache line contains 4 doubles as sizeof(double) is 8 bytes
#elif defined MGONGPU_FPTYPE_FLOAT
  const int neppM = 8; // DEFAULT: one 32-byte cache line contains 8 floats as sizeof(float) is 4 bytes
#endif
  //const int neppM = 1;  // *** NB: this is equivalent to AOS ***
  //const int neppM = 32; // older default

  // Number of Events Per Page in the random number AOSOA memory layout
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
  const int neppR = 8; // HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
  //const int neppR = 1; // AOS (tests of sectors/requests)

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
// Arguments (not used so far): text is __FUNCTION__, code is 0 (start) or 1 (end)
#if defined __CUDACC__ && defined MGONGPU_NSIGHT_DEBUG
#define mgDebugDeclare() \
  __shared__ float mgDebugCounter[mgOnGpu::ntpbMAX];
#define mgDebugInitialise() \
  { mgDebugCounter[threadIdx.x]=0; }
#define mgDebug( code, text ) \
  { mgDebugCounter[threadIdx.x] += 1; }
#define mgDebugFinalise() \
  { if ( blockIdx.x == 0 && threadIdx.x == 0 ) printf( "MGDEBUG: counter=%%f\n", mgDebugCounter[threadIdx.x] ); }
#else
#define mgDebugDeclare() \
  /*noop*/
#define mgDebugInitialise() \
  { /*noop*/ }
#define mgDebug( code, text ) \
  { /*noop*/ }
#define mgDebugFinalise() \
  { /*noop*/ }
#endif

// Define empty CUDA declaration specifiers for C++
#ifndef __CUDACC__
#define __global__
//#define __host__
#define __device__
#endif

#endif // MGONGPUCONFIG_H
