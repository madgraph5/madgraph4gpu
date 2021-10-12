//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

//#include <cmath>
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm
{

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
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ipzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void imzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void vxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv* vc,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void sxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype,                   // WARNING: "smass" unused (missing in Fortran)
               const int,                      // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[3],                // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void opzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void omzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1P0_3( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;    // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_3( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;      // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_3( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;      // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_0( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype_sv V3[],              // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* vertex ) ALWAYS_INLINE; // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_3( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;    // output: wavefunction3[6]

  //--------------------------------------------------------------------------

} // end namespace MG5_sm

#endif // HelAmps_sm_H

