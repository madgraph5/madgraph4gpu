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

  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void ixxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void ixxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void ipzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void ipzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void imzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void imzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void ixzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__ INLINE
  void ixzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void vxxxxx( const fptype_sv* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void sxxxxx( const fptype_sv* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void opzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void opzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void omzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS>
  __host__ __device__ INLINE
  void oxzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ INLINE
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

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

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ INLINE
  void FFV2_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ INLINE
  void FFV2_3( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ INLINE
  void FFV4_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ INLINE
  void FFV4_3( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ INLINE
  void FFV2_4_0( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype_sv V3[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ INLINE
  void FFV2_4_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;

  //==========================================================================

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS>
  __host__ __device__
  void ixxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    if ( fmass != 0. )
    {
      const fptype_sv pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
#ifndef MGONGPU_CPPSIMD
      if ( pp == 0. )
      {
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
        //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const int ip = ( 1 + nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        const int im = ( 1 - nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        fi[2] = cxmake( ip * sqm[ip], 0 );
        fi[3] = cxmake( im * nsf * sqm[ip], 0 );
        fi[4] = cxmake( ip * nsf * sqm[im], 0 );
        fi[5] = cxmake( im * sqm[im], 0 );
      }
      else
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype omega[2] = { fpsqrt( pvec0 + pp ), 0. };
        omega[1] = fmass / omega[0];
        const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype pp3 = fpmax( pp + pvec3, 0. );
        const cxtype chi[2] = { cxmake( fpsqrt ( pp3 * (fptype)0.5 / pp ), 0. ),
                                ( pp3 == 0. ?
                                  cxmake( -nh, 0. ) :
                                  cxmake( nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        fi[2] = sfomega[0] * chi[im];
        fi[3] = sfomega[0] * chi[ip];
        fi[4] = sfomega[1] * chi[im];
        fi[5] = sfomega[1] * chi[ip];
      }
#else
      const int ip = ( 1 + nh ) / 2;
      const int im = ( 1 - nh ) / 2;
      // Branch A: pp == 0.
      // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
      fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0 }; // possibility of negative fermion masses (NB: SCALAR!)
      sqm[1] = ( fmass < 0 ? -sqm[0] : sqm[0] ); // AV: removed an abs here (as above)
      const cxtype fiA_2 = ip * sqm[ip]; // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_3 = im * nsf * sqm[ip]; // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_4 = ip * nsf * sqm[im]; // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_5 = im * sqm[im]; // scalar cxtype: real part initialised from fptype, imag part = 0
      // Branch B: pp != 0.
      const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                             fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
      fptype_v omega[2] = { fpsqrt( pvec0 + pp ), 0 };
      omega[1] = fmass / omega[0];
      const fptype_v sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const fptype_v pp3 = fpmax( pp + pvec3, 0 );
      const cxtype_v chi[2] = { cxmake( fpsqrt ( pp3 * 0.5 / pp ), 0 ),
                                cxternary( ( pp3 == 0. ),
                                           cxmake( -nh, 0 ),
                                           cxmake( (fptype)nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
      const cxtype_v fiB_2 = sfomega[0] * chi[im];
      const cxtype_v fiB_3 = sfomega[0] * chi[ip];
      const cxtype_v fiB_4 = sfomega[1] * chi[im];
      const cxtype_v fiB_5 = sfomega[1] * chi[ip];
      // Choose between the results from branch A and branch B
      const bool_v mask = ( pp == 0. );
      fi[2] = cxternary( mask, fiA_2, fiB_2 );
      fi[3] = cxternary( mask, fiA_3, fiB_3 );
      fi[4] = cxternary( mask, fiA_4, fiB_4 );
      fi[5] = cxternary( mask, fiA_5, fiB_5 );
#endif
    }
    else
    {
      const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                          fptype_sv{0}, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
      const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ), cxternary( ( sqp0p3 == 0. ),
                                                                  cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0 ), 0. ),
                                                                  cxmake( (fptype)nh * pvec1, pvec2 ) / sqp0p3 ) };
      if ( nh == 1 )
      {
        fi[2] = cxzero_sv();
        fi[3] = cxzero_sv();
        fi[4] = chi[0];
        fi[5] = chi[1];
      }
      else
      {
        fi[2] = chi[1];
        fi[3] = chi[0];
        fi[4] = cxzero_sv();
        fi[5] = cxzero_sv();
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS>
  __host__ __device__
  void ipzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    fi[0] = cxmake( -pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv sqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
    fi[2] = fi[1];
    if( nh == 1 )
    {
      fi[3] = fi[1];
      fi[4] = sqp0p3;
    }
    else
    {
      fi[3] = sqp0p3;
      fi[4] = fi[1];
    }
    fi[5] = fi[1];
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  template<class M_ACCESS>
  __host__ __device__
  void imzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    fi[0] = cxmake( pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv chi = cxmake( -(fptype)nhel * fpsqrt( -2. * pvec3 ), 0. );
    fi[3] = cxzero_sv();
    fi[4] = cxzero_sv();
    if ( nh == 1 )
    {
      fi[2] = cxzero_sv();
      fi[5] = chi;
    }
    else
    {
      fi[2] = chi;
      fi[5] = cxzero_sv();
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  template<class M_ACCESS>
  __host__ __device__
  void ixzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    //fi[0] = cxmake( -pvec0 * nsf, -pvec2 * nsf ); // AV: BUG! not the same as ixxxxx
    //fi[1] = cxmake( -pvec0 * nsf, -pvec1 * nsf ); // AV: BUG! not the same as ixxxxx
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf ); // AV: BUG FIX
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf ); // AV: BUG FIX
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1/sqp0p3, pvec2/sqp0p3 );
    if ( nh == 1 )
    {
      fi[2] = cxzero_sv();
      fi[3] = cxzero_sv();
      fi[4] = chi0;
      fi[5] = chi1;
    }
    else
    {
      fi[2] = chi1;
      fi[3] = chi0;
      fi[4] = cxzero_sv();
      fi[5] = cxzero_sv();
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS>
  __host__ __device__
  void opzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    fo[0] = cxmake( pvec3 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv csqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
    fo[3] = cxzero_sv();
    fo[4] = cxzero_sv();
    if ( nh == 1 )
    {
      fo[2] = csqp0p3;
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[5] = csqp0p3;
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  template<class M_ACCESS>
  __host__ __device__
  void oxzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
    if ( nh == 1 )
    {
      fo[2] = chi0;
      fo[3] = chi1;
      fo[4] = cxzero_sv();
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[3] = cxzero_sv();
      fo[4] = chi1;
      fo[5] = chi0;
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace MG5_sm

#endif // HelAmps_sm_H

