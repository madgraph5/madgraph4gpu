//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm
{

  // =============================================================================
  // *** Generic pattern: kernelAccessFunction( buffer, additional_indexes ) ***
  // =============================================================================

  // Kernel access function (WITHOUT an explicit event number) for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the 4-momenta for one event or one SIMD vector of events
  // (Non-const memory access)
  __device__ inline
  fptype_sv& kernelAccessMomenta( fptype_sv* buffer,
                                  const int ip4
#ifdef __CUDACC__
                                  , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                  );

  // (Const memory access)
  __device__ inline
  const fptype_sv& kernelAccessConstMomenta( const fptype_sv* buffer,
                                             const int ip4
#ifdef __CUDACC__
                                             , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                             );

  // Kernel access function (WITHOUT an explicit event number) for amplitudes
  // Input: a memory buffer for an arbitrary number of events
  // Output: the amplitude for one event or one SIMD vector of events
  // (Non-const memory access)
  __device__ inline
  cxtype_sv& kernelAccessAmplitudes( cxtype_sv* buffer );

  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ipzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void imzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void vxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv* vc                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void sxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype,                   // WARNING: "smass" unused (missing in Fortran)
               const int,                      // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[3]                 // output: wavefunction[3] - not [6], this is for scalars
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void opzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void omzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo                   // output: wavefunction[(nw6==6)]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1_0( const cxtype_sv F1[],    // input: wavefunction1[6]
               const cxtype_sv F2[],    // input: wavefunction2[6]
               const cxtype_sv V3[],    // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* amplitudes )  // output: amplitude
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1P0_3( const cxtype_sv F1[],  // input: wavefunction1[6]
                 const cxtype_sv F2[],  // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )       // output: wavefunction3[6]
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_0( const cxtype_sv F1[],    // input: wavefunction1[6]
               const cxtype_sv F2[],    // input: wavefunction2[6]
               const cxtype_sv V3[],    // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* amplitudes )  // output: amplitude
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_3( const cxtype_sv F1[],    // input: wavefunction1[6]
               const cxtype_sv F2[],    // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] )         // output: wavefunction3[6]
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_0( const cxtype_sv F1[],     // input: wavefunction1[6]
               const cxtype_sv F2[],     // input: wavefunction2[6]
               const cxtype_sv V3[],     // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* amplitudes )   // output: amplitude
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_3( const cxtype_sv F1[],     // input: wavefunction1[6]
               const cxtype_sv F2[],     // input: wavefunction2[6]
               const cxtype_sv COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] )          // output: wavefunction3[6]
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_0( const cxtype_sv F1[],    // input: wavefunction1[6]
                 const cxtype_sv F2[],    // input: wavefunction2[6]
                 const cxtype_sv V3[],    // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* amplitudes )  // output: amplitude
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_3( const cxtype_sv F1[],    // input: wavefunction1[6]
                 const cxtype_sv F2[],    // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )         // output: wavefunction3[6]
    ALWAYS_INLINE;

  //--------------------------------------------------------------------------

}

#endif // HelAmps_sm_H
