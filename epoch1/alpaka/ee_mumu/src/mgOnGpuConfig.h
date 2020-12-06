#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// ** NB1 Throughputs (e.g. 6.8E8) are events/sec for "./gcheck.exe -p 65536 128 12"
// ** NB2 Baseline on b7g47n0004 fluctuates (probably depends on load on other VMs)

// Memory layout for momenta: AOSOA, AOS, SOA (CHOOSE ONLY ONE)
// AOSOA (ASA) layout is hardcoded: fine-tune it using the nepopR and neppM parameters below

// Random generator type: This is curand on the standard cupla source. alsimple is a
// small untility class intended for use with the Alpaka port.
//#define MGONGPU_RANDTYPE_CURAND 1 // original default
#define MGONGPU_RANDTYPE_ALSIMPLE 1

// Random number generation location (CHOOSE ONLY ONE)
//#define MGONGPU_RAND_ONDEVICE 1 // original default, but currently not available for alsimple
#define MGONGPU_RAND_ONHOST 1

// Memory choice for wavefunctions: registries/"local", global, shared (CHOOSE ONLY ONE)
// Local storage (registries plus spillover to local) is hardcoded: fine tune it using maxrregcount in the Makefile
// [NB: new throughputs on 1GPU/4CPU system]

// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default (~6.8E8)
//#define MGONGPU_FPTYPE_FLOAT 1 // 2.4x faster (~1.64E9 against 6.8E8)

// Complex type: thrust, cucomplex or alsimple (CHOOSE ONLY ONE)
//#define MGONGPU_CXTYPE_THRUST 1 // original default (~6.8E8)
//#define MGONGPU_CXTYPE_CUCOMPLEX 1 // ~5% slower (6.5E8 against 6.8E8)
#define MGONGPU_CXTYPE_ALSIMPLE 1

// Cuda nsight compute (ncu) debug: add dummy lines to ease SASS program flow navigation
#undef MGONGPU_NSIGHT_DEBUG // default
//#define MGONGPU_NSIGHT_DEBUG 1

namespace mgOnGpu
{
  // --- Physics process-specific constants that are best declared at compile time

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
  const int nw6 = 6; // dimension of each wavefunction (see KEK 91-11)

  const int npari = 2; // #particles in the initial state (incoming): e+ e-
  const int nparf = 2; // #particles in the final state (outgoing): mu+ mu-
  const int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  const int nwf = 5; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)

  const int ncomb = 16; // #helicity combinations: 16=2(spin up/down for fermions)**4(npar)

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  const int ntpbMAX = 256;

  // Number of Events Per Page in the random number AOSOA (ASA) structure
  // (this is best kept as a compile-time constant: see issue #23)
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
#if defined MGONGPU_FPTYPE_DOUBLE
//  const int neppR = 4; // DEFAULT: one 32-byte cache line contains 4 doubles as sizeof(double) is 8 bytes
#elif defined MGONGPU_FPTYPE_FLOAT
//  const int neppR = 8; // DEFAULT: one 32-byte cache line contains 8 floats as sizeof(float) is 4 bytes
#endif
const int neppR = 1;  // *** NB: this is equivalent to AOS ***
  //const int neppR = 32; // older default

  // Number of Events Per Page in the momenta AOSOA (ASA) structure
  // (this is best kept as a compile-time constant: see issue #23)
#if defined MGONGPU_FPTYPE_DOUBLE
//  const int neppM = 4; // DEFAULT: one 32-byte cache line contains 4 doubles as sizeof(double) is 8 bytes
#elif defined MGONGPU_FPTYPE_FLOAT
//  const int neppM = 8; // DEFAULT: one 32-byte cache line contains 8 floats as sizeof(float) is 4 bytes
#endif
const int neppM = 1;  // *** NB: this is equivalent to AOS ***
  //const int neppM = 32; // older default

}

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
  { if ( blockIdx.x == 0 && threadIdx.x == 0 ) printf( "MGDEBUG: counter=%f\n", mgDebugCounter[threadIdx.x] ); }
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

#endif // MGONGPUCONFIG_H
