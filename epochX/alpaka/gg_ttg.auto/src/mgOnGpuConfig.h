#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// ** NB1 Throughputs (e.g. 6.8E8) are events/sec for "./gcheck.exe -p 65536 128 12"
// ** NB2 Baseline on b7g47n0004 fluctuates (probably depends on load on other VMs)

// Choose how random numbers are generated
// If one of these macros has been set from outside with e.g. -DMGONGPU_COMMONRAND_ONHOST, nothing happens
#if not defined MGONGPU_CURAND_ONDEVICE and not defined MGONGPU_CURAND_ONHOST and not defined MGONGPU_COMMONRAND_ONHOST
#ifdef ALPAKA
#define MGONGPU_COMMONRAND_ONHOST 1 // default for alpaka (common rand: on host with alpaka)
#else
// Curand random number generation (CHOOSE ONLY ONE)
#define MGONGPU_CURAND_ONDEVICE 1 // default for non-alpaka (curand: CUDA on device, C++ on host)
//#define MGONGPU_CURAND_ONHOST 1 // (curand: CUDA on host, C++ on host)
//#define MGONGPU_COMMONRAND_ONHOST 1 // (common rand: CUDA on host, C++ on host)
#endif
#endif

// Choose floating point precision
// If one of these macros has been set from outside with e.g. -DMGONGPU_FPTYPE_FLOAT, nothing happens (issue #167)
#if not defined MGONGPU_FPTYPE_DOUBLE and not defined MGONGPU_FPTYPE_FLOAT
// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default (~6.8E8)
//#define MGONGPU_FPTYPE_FLOAT 1 // 2.4x faster (~1.64E9 against 6.8E8)
#endif

// Choose whether to inline all HelAmps functions
// This optimization can gain almost a factor 4 in C++, similar to -flto (issue #229)
// By default, do not inline, but allow this macros to be set from outside with e.g. -DMGONGPU_INLINE_HELAMPS
//#undef MGONGPU_INLINE_HELAMPS // default
////#define MGONGPU_INLINE_HELAMPS 1

// Complex type in cuda: thrust or cucomplex (CHOOSE ONLY ONE)
#ifdef ALPAKA
// for alpaka by default use custom complex class that doesn't depend on cuda or on thrust
#define MGONGPU_CXTYPE_ALSIMPLE 1
#elif defined(__CUDACC__)
#define MGONGPU_CXTYPE_THRUST 1 // default (~6.8E8)
//#define MGONGPU_CXTYPE_CUCOMPLEX 1 // ~5 percent slower (6.5E8 against 6.8E8)
#endif

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
#ifdef __CUDACC__
#undef MGONGPU_NSIGHT_DEBUG // default
//#define MGONGPU_NSIGHT_DEBUG 1
#endif

// SANITY CHECKS (random numbers)
#if defined MGONGPU_CURAND_ONDEVICE and defined MGONGPU_CURAND_ONHOST
#error You must CHOOSE ONLY ONE of MGONGPU_CURAND_ONDEVICE, MGONGPU_CURAND_ONHOST or MGONGPU_COMMONRAND_ONHOST
#elif defined MGONGPU_CURAND_ONDEVICE and defined MGONGPU_COMMONRAND_ONHOST
#error You must CHOOSE ONLY ONE of MGONGPU_CURAND_ONDEVICE, MGONGPU_CURAND_ONHOST or MGONGPU_COMMONRAND_ONHOST
#elif defined MGONGPU_CURAND_ONHOST and defined MGONGPU_COMMONRAND_ONHOST
#error You must CHOOSE ONLY ONE of MGONGPU_CURAND_ONDEVICE, MGONGPU_CURAND_ONHOST or MGONGPU_COMMONRAND_ONHOST
#endif

// SANITY CHECKS (floating point precision)
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE_FLOAT
#error You must CHOOSE ONLY ONE of MGONGPU_FPTYPE_DOUBLE or defined MGONGPU_FPTYPE_FLOAT
#endif

// SANITY CHECKS (complex number implementation)
#if defined(ALPAKA) || defined(__CUDACC__)
#if defined MGONGPU_CXTYPE_THRUST and defined MGONGPU_CXTYPE_CUCOMPLEX
#error You must CHOOSE ONLY ONE of MGONGPU_CXTYPE_THRUST, MGONGPU_CXTYPE_CUCOMPLEX or MGONGPU_CXTYPE_ALSIMPLE
#endif
#if defined MGONGPU_CXTYPE_THRUST and defined MGONGPU_CXTYPE_ALSIMPLE
#error You must CHOOSE ONLY ONE of MGONGPU_CXTYPE_THRUST, MGONGPU_CXTYPE_CUCOMPLEX or MGONGPU_CXTYPE_ALSIMPLE
#endif
#if defined MGONGPU_CXTYPE_CUCOMPLEX and defined MGONGPU_CXTYPE_ALSIMPLE
#error You must CHOOSE ONLY ONE of MGONGPU_CXTYPE_THRUST, MGONGPU_CXTYPE_CUCOMPLEX or MGONGPU_CXTYPE_ALSIMPLE
#endif
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

  const int np4 = 4; // dimensions of 4-momenta (E,px,py,pz)

  const int npari = 2; // #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
  const int nparf = 3; // #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
  const int npar = npari + nparf; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-

  const int ncomb = 32; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  const int nw6 = 6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  const int nwf = 12; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  //const int ntpbMAX = 256; // AV Apr2021: why had I set this to 256?
  const int ntpbMAX = 1024; // NB: 512 is ok, but 1024 does fail with "too many resources requested for launch"

  // Vector sizes for AOSOA memory layouts (GPU coalesced memory access, CPU SIMD vectorization)
  // (these are all best kept as a compile-time constants: see issue #23)

  // Number of Events Per Page in the momenta AOSOA memory layout
#if defined(ALPAKA) || defined(__CUDACC__)
#undef MGONGPU_CPPSIMD
  // -----------------------------------------------------------------------------------------
  // --- GPUs: neppM must be a power of 2 times the number of fptype's in a 32-byte cacheline
  // --- This is relevant to ensure coalesced access to momenta in global memory
  // --- Note that neppR is hardcoded and may differ from neppM and neppV on some platforms
  // -----------------------------------------------------------------------------------------
  // AOS needed for some Alpaka configurations
#ifdef ALPAKA_AOS
  const int neppM = 1;
#else
  //const int neppM = 64/sizeof(fptype); // 2x 32-byte GPU cache lines: 8 (DOUBLE) or 16 (FLOAT)
  const int neppM = 32/sizeof(fptype); // (DEFAULT) 32-byte GPU cache line: 4 (DOUBLE) or 8 (FLOAT)
  //const int neppM = 1;  // *** NB: this is equivalent to AOS ***
  //const int neppM = 32; // older default
#endif
#else
  // -----------------------------------------------------------------------------------------
  // --- CPUs: neppM must be exactly equal to the number of fptype's in a vector register
  // --- [DEFAULT is 256-width AVX512 aka "512y" on CPUs, 32-byte as GPUs, faster than AVX2]
  // --- The logic of the code requires the size neppV of fptype_v to be equal to neppM
  // --- Note that neppR is hardcoded and may differ from neppM and neppV on some platforms
  // -----------------------------------------------------------------------------------------
#if defined __AVX512VL__
#define MGONGPU_CPPSIMD 1
#ifdef MGONGPU_PVW512
  const int neppM = 64/sizeof(fptype); // "512z" AVX512 with 512 width (512-bit ie 64-byte): 8 (DOUBLE) or 16 (FLOAT)
#else
  const int neppM = 32/sizeof(fptype); // "512y" AVX512 with 256 width (256-bit ie 32-byte): 4 (DOUBLE) or 8 (FLOAT) [gcc DEFAULT]
#endif
#elif defined __AVX2__
#define MGONGPU_CPPSIMD 1
  const int neppM = 32/sizeof(fptype); // "avx2" AVX2 (256-bit ie 32-byte): 4 (DOUBLE) or 8 (FLOAT) [clang DEFAULT]
#elif defined __SSE4_2__
#define MGONGPU_CPPSIMD 1
  const int neppM = 16/sizeof(fptype); // "sse4" SSE4.2 (128-bit ie 16-byte): 2 (DOUBLE) or 4 (FLOAT)
#else
#undef MGONGPU_CPPSIMD
  const int neppM = 1;  // "none" i.e. no SIMD (*** NB: this is equivalent to AOS ***)
#endif
#endif

  // Number of Events Per Page in the random number AOSOA memory layout
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
  // AOS needed for some Alpaka configurations
#ifdef ALPAKA_AOS
  const int neppR = 1;
#else
  const int neppR = 8; // HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
  //const int neppR = 1; // AOS (tests of sectors/requests)
#endif

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
// Arguments (not used so far): text is __FUNCTION__, code is 0 (start) or 1 (end)
#if defined __CUDACC__ && defined MGONGPU_NSIGHT_DEBUG
#define mgDebugDeclare()                              \
  __shared__ float mgDebugCounter[mgOnGpu::ntpbMAX];
#define mgDebugInitialise()                     \
  { mgDebugCounter[threadIdx.x]=0; }
#define mgDebug( code, text )                   \
  { mgDebugCounter[threadIdx.x] += 1; }
#define mgDebugFinalise()                                               \
  { if ( blockIdx.x == 0 && threadIdx.x == 0 ) printf( "MGDEBUG: counter=%f\n", mgDebugCounter[threadIdx.x] ); }
#else
#define mgDebugDeclare()                        \
  /*noop*/
#define mgDebugInitialise()                     \
  { /*noop*/ }
#define mgDebug( code, text )                   \
  { /*noop*/ }
#define mgDebugFinalise()                       \
  { /*noop*/ }
#endif

// Define empty CUDA declaration specifiers for C++
#if !defined(ALPAKA) && !defined(__CUDACC__)
#define __global__
//#define __host__
#define __device__
#endif

#endif // MGONGPUCONFIG_H
