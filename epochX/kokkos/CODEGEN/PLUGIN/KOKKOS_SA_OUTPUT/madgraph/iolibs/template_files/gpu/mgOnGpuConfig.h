#ifndef MGONGPUCONFIG_H
#define MGONGPUCONFIG_H 1

// Choose floating point precision
// If one of these macros has been set from outside with e.g. -DMGONGPU_FPTYPE_FLOAT, nothing happens (issue #167)
#if not defined MGONGPU_FPTYPE_DOUBLE and not defined MGONGPU_FPTYPE_FLOAT
// Floating point precision (CHOOSE ONLY ONE)
#define MGONGPU_FPTYPE_DOUBLE 1 // default (~6.8E8)
//#define MGONGPU_FPTYPE_FLOAT 1 // 2.4x faster (~1.64E9 against 6.8E8)
#endif

// Choose whether to inline all HelAmps functions
// This optimization can gain almost a factor 4 in C++, similar to -flto (issue #229)
// By default, do not inline, but allow this macro to be set from outside with e.g. -DMGONGPU_INLINE_HELAMPS
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

  const int np4 = 4; // dimensions of 4-momenta (E,px,py,pz)

  const int npari = %(nincoming)d; // #particles in the initial state (incoming): e.g. 2 (e+ e-) for e+ e- -> mu+ mu-
  const int nparf = %(noutcoming)d; // #particles in the final state (outgoing): e.g. 2 (mu+ mu-) for e+ e- -> mu+ mu-
  const int npar = npari + nparf; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-

  const int ncomb = %(nbhel)d; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  const int nw6 = %(wavefuncsize)d; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  const int nwf = %(nwavefunc)d; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  
  const int ncouplings = %(ncouplings)d;
  const int ncouplingstimes2 = %(ncouplingstimes2)d;
  const int nparams = %(nparams)d;

  // --- Platform-specific software implementation details

  // Maximum number of blocks per grid
  // ** NB Some arrays of pointers will be allocated statically to fit all these blocks
  // ** (the actual memory for each block will then be allocated dynamically only for existing blocks)
  //const int nbpgMAX = 2048;

  // Maximum number of threads per block
  //const int ntpbMAX = 256; // AV Apr2021: why had I set this to 256?
#ifdef MGONGPU_NTPBMAX
  const int ntpbMAX = MGONGPU_NTPBMAX;
#else
  const int ntpbMAX = 256; // NB: 512 is ok, but 1024 does fail with "too many resources requested for launch"
#endif
// -----------------------------------------------------------------------------------------------
  // --- GPUs: neppM is best set to a power of 2 times the number of fptype's in a 32-byte cacheline
  // --- This is relevant to ensure coalesced access to momenta in global memory
  // --- Note that neppR is hardcoded and may differ from neppM and neppV on some platforms
  // -----------------------------------------------------------------------------------------------
#ifdef MGONGPU_NEPPM
  static constexpr int  neppM = MGONGPU_NEPPM;
#else
  //static constexpr int neppM = 64/sizeof(fptype); // 2x 32-byte GPU cache lines (512 bits): 8 (DOUBLE) or 16 (FLOAT)
  static constexpr int neppM = 32/sizeof(fptype); // (DEFAULT) 32-byte GPU cache line (256 bits): 4 (DOUBLE) or 8 (FLOAT)
  //static constexpr int neppM = 1;  // *** NB: this is equivalent to AOS *** (slower: 1.03E9 instead of 1.11E9 in eemumu)
#endif
}
// Expose typedefs and operators outside the namespace
using mgOnGpu::fptype;

#endif // MGONGPUCONFIG_H



/* #ifndef MGKOKKOSCONFIG_H */
/* #define MGKOKKOSCONFIG_H 1 */

/* namespace mgKokkos */
/* { */
/*   // --- Physics process-specific constants that are best declared at compile time */

/*   constexpr int np4 = 4; // the dimension of 4-momenta (E,px,py,pz) */
/*   constexpr int nw6 = %(wavefuncsize)d; // dimension of each wavefunction (see KEK 91-11) */

/*   constexpr int npari = %(nincoming)d; // #particles in the initial state (incoming): e+ e- */
/*   constexpr int nparf = %(noutcoming)d; // #particles in the final state (outgoing): mu+ mu- */
/*   constexpr int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu- */

/*   constexpr int nwf = %(nwavefunc)d; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z) */

/*   constexpr int ncomb = %(nbhel)d; // #helicity combinations: 16=2(spin up/down for fermions)**4(npar) */

/*   // --- Platform-specific software implementation details */

/*   // Maximum number of threads per block */
/*   constexpr int ntpbMAX = 256; */

/* } */

/* #endif // MGKOKKOSCONFIG_H */
