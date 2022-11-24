#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

#define MGONGPU_NDCOUP 3
#define MGONGPU_NICOUP 0

// HARDCODED AT CODE GENERATION TIME: DO NOT MODIFY (#473)
// There are two different code bases for standalone_sycl (without multichannel) and madevent+sycl (with multichannel)
#define MGONGPU_SUPPORTS_MULTICHANNEL 1

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
#ifdef __CUDACC__
#define MGONGPU_CXTYPE_THRUST 1 // default (~6.8E8)
//#define MGONGPU_CXTYPE_CUCOMPLEX 1 // ~5 percent slower (6.5E8 against 6.8E8)
#endif

// SANITY CHECKS (floating point precision)
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE_FLOAT
#error You must CHOOSE ONLY ONE of MGONGPU_FPTYPE_DOUBLE or defined MGONGPU_FPTYPE_FLOAT
#endif

// SANITY CHECKS (complex number implementation)
#ifdef __CUDACC__
#if defined MGONGPU_CXTYPE_THRUST and defined MGONGPU_CXTYPE_CUCOMPLEX
#error You must CHOOSE ONLY ONE of MGONGPU_CXTYPE_THRUST or MGONGPU_CXTYPE_CUCOMPLEX
#endif
#endif

namespace mgOnGpu
{

  // --- Type definitions

  // Floating point type: fptype
#if defined MGONGPU_FPTYPE_DOUBLE
  typedef double fptype; // double precision (8 bytes, fp64)
#define FPTYPE_NAME "DOUBLE"
#elif defined MGONGPU_FPTYPE_FLOAT
  typedef float fptype; // single precision (4 bytes, fp32)
#define FPTYPE_NAME "FLOAT"
#endif

  // --- Physics process-specific constants that are best declared at compile time

  static constexpr int np4 = 4; // dimensions of 4-momenta (E,px,py,pz)

  static constexpr int npari = 2; // #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
  static constexpr int nparf = 3; // #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
  static constexpr int npar = npari + nparf; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-

  static constexpr int ncomb = 32; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  static constexpr int nw6 = 6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  static constexpr int nwf = 12; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  
  static constexpr int ncouplings = 3;
  static constexpr int ncouplingstimes2 = 6;
  static constexpr int nparams = 2;

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  //const int ntpbMAX = 256; // AV Apr2021: why had I set this to 256?
#ifdef MGONGPU_NTPBMAX
  static constexpr int ntpbMAX = MGONGPU_NTPBMAX;
#else
  static constexpr int ntpbMAX = 1024; // NB: 512 is ok, but 1024 does fail with "too many resources requested for launch"
#endif

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;


#endif // MGONGPUCONFIG_H
