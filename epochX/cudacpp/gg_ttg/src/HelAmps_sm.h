//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
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

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void ixxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void ipzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void imzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__ INLINE
  void ixzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void vxxxxx( const fptype_sv* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void sxxxxx( const fptype_sv* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void opzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void omzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
  __device__ INLINE
  void VVV1_0( const cxtype_sv V1[],
               const cxtype_sv V2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  __device__ INLINE
  void VVV1P0_1( const cxtype_sv V2[],
                 const cxtype_sv V3[],
                 const cxtype COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ INLINE
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  __device__ INLINE
  void FFV1_1( const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M1,
               const fptype W1,
               cxtype_sv F1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  __device__ INLINE
  void FFV1_2( const cxtype_sv F1[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M2,
               const fptype W2,
               cxtype_sv F2[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ INLINE
  void FFV1P0_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  __device__ INLINE
  void VVVV1P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  __device__ INLINE
  void VVVV3P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  __device__ INLINE
  void VVVV4P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

} // end namespace MG5_sm

#endif // HelAmps_sm_H

