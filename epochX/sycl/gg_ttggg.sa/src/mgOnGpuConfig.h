#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

#define MGONGPU_NDCOUP 3
#define MGONGPU_NICOUP 0
#define MGONGPU_FOURVECTOR_DIM 4
#define NPAR 7


//Sets vectorization level when using vectorizable complex types
#if not defined MGONGPU_VEC_DIM
    #define MGONGPU_VEC_DIM 1
#endif
#if \
    MGONGPU_VEC_DIM !=  1 || \
    MGONGPU_VEC_DIM !=  2 || \
    MGONGPU_VEC_DIM !=  3 || \
    MGONGPU_VEC_DIM !=  4 || \
    MGONGPU_VEC_DIM !=  8 || \
    MGONGPU_VEC_DIM != 16

    #error You must set MGONGPU_VEC_DIM to 1, 2, 3, 4, 8, or 16.
#endif

// Set complex number library
#if\
    not defined MGONGPU_COMPLEX_CXSMPL    || \
    not defined MGONGPU_COMPLEX_EXTRAS    || \
    not defined MGONGPU_COMPLEX_STD       || \
    not defined MGONGPU_COMPLEX_ONEAPI    || \
    not defined MGONGPU_COMPLEX_CUTHRUST  || \
    not defined MGONGPU_COMPLEX_CUCOMPLEX

    #define MGONGPU_COMPLEX_STD 1
#endif

#if\
    defined() + \
    defined(MGONGPU_COMPLEX_CXSMPL)   + \
    defined(MGONGPU_COMPLEX_EXTRAS)   + \
    defined(MGONGPU_COMPLEX_STD)      + \
    defined(MGONGPU_COMPLEX_ONEAPI)   + \
    defined(MGONGPU_COMPLEX_CUTHRUST) + \
    defined(MGONGPU_COMPLEX_CUCOMPLEX)  \
    != 1

    #error You must CHOOSE (ONE AND) ONLY ONE complex number library
#endif

// HARDCODED AT CODE GENERATION TIME: DO NOT MODIFY (#473)
// There are two different code bases for standalone_sycl (without multichannel) and madevent+sycl (with multichannel)
#undef MGONGPU_SUPPORTS_MULTICHANNEL

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

// SANITY CHECKS (floating point precision)
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE_FLOAT
#error You must CHOOSE (ONE AND) ONLY ONE of MGONGPU_FPTYPE_DOUBLE or defined MGONGPU_FPTYPE_FLOAT
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

  static constexpr unsigned int np4 = 4; // dimensions of 4-momenta (E,px,py,pz)

  static constexpr unsigned int npari = 2; // #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
  static constexpr unsigned int nparf = 5; // #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
  static constexpr unsigned int npar = npari + nparf; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-

  static constexpr unsigned int ncomb = 128; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  static constexpr unsigned int nw6 = 6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  static constexpr unsigned int nwf = 121; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  
  static constexpr unsigned int ncouplings = 3;
  static constexpr unsigned int ncouplingstimes2 = 6;
  static constexpr unsigned int nparams = 2;

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  //const int ntpbMAX = 256; // AV Apr2021: why had I set this to 256?
#ifdef MGONGPU_NTPBMAX
  static constexpr unsigned int ntpbMAX = MGONGPU_NTPBMAX;
#else
  static constexpr unsigned int ntpbMAX = 1024; // NB: 512 is ok, but 1024 does fail with "too many resources requested for launch"
#endif

  // -----------------------------------------------------------------------------------------------
  // --- GPUs: neppM is best set to a power of 2 times the number of fptype's in a 32-byte cacheline
  // --- This is relevant to ensure coalesced access to momenta in global memory
  // --- Note that neppR is hardcoded and may differ from neppM and neppV on some platforms
  // -----------------------------------------------------------------------------------------------
  static constexpr unsigned int neppM = 1; // FIXME neppM unused now, needs to be removed from random number generation

  // Number of Events Per Page in the random number AOSOA memory layout
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
  static constexpr unsigned int neppR = 8; // HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
  //const int neppR = 1; // AOS (tests of sectors/requests)

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;

// SYCL debug: add dummy lines to ease SASS program flow navigation
// Arguments (not used so far): text is __FUNCTION__, code is 0 (start) or 1 (end)
#if defined SYCL_DEBUG
#define mgDebugDeclare()                              \
  //__shared__ float mgDebugCounter[mgOnGpu::ntpbMAX];
#define mgDebugInitialise()                     \
  { /*mgDebugCounter[threadIdx.x] = 0;*/ }
#define mgDebug( code, text )                   \
  { /*mgDebugCounter[threadIdx.x] += 1; */}
#define mgDebugFinalise()                                               \
  { /*if ( blockIdx.x == 0 && threadIdx.x == 0 ) printf( "MGDEBUG: counter=%f\n", mgDebugCounter[threadIdx.x] ); */}
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


#endif // MGONGPUCONFIG_H
