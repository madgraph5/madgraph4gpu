//==========================================================================
// This file has been automatically generated for SYCL/C++ standalone by
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

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
   INLINE
  void VVV1_0( const cxtype_sv V1[],
               const cxtype_sv V2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
   INLINE
  void VVV1P0_1( const cxtype_sv V2[],
                 const cxtype_sv V3[],
                 const cxtype COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
   INLINE
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
   INLINE
  void FFV1_1( const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M1,
               const fptype W1,
               cxtype_sv F1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
   INLINE
  void FFV1_2( const cxtype_sv F1[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M2,
               const fptype W2,
               cxtype_sv F2[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
   INLINE
  void FFV1P0_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
   INLINE
  void VVVV1_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
   INLINE
  void VVVV1P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
   INLINE
  void VVVV3_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
   INLINE
  void VVVV3P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
   INLINE
  void VVVV4_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
   INLINE
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

