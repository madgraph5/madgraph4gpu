// Copyright (C) 2010 The ALOHA Development team and Contributors.
// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Sep 2010) for the MG5aMC backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.5.2_lo_vect, 2023-11-08
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuVectors.h"

#include "Parameters_sm.h"

#include <cassert>
//#include <cmath>
//#include <cstdlib>
//#include <iomanip>
//#include <iostream>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__( ( always_inline ) )
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  ixxxxx( const fptype momenta[], // input: momenta
          const fptype fmass,     // input: fermion mass
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  ipzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  imzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  ixzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  vxxxxx( const fptype momenta[], // input: momenta
          const fptype vmass,     // input: vector boson mass
          const int nhel,         // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
          const int nsv,          // input: +1 (final) or -1 (initial)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  sxxxxx( const fptype momenta[], // input: momenta
          //const fptype,                 // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
          //const int,                    // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
          const int nss,          // input: +1 (final) or -1 (initial)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  oxxxxx( const fptype momenta[], // input: momenta
          const fptype fmass,     // input: fermion mass
          const int nhel,         // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  opzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  omzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ INLINE void
  oxzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar          // input: particle# out of npar
          ) ALWAYS_INLINE;

  //==========================================================================

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  ixxxxx( const fptype momenta[], // input: momenta
          const fptype fmass,     // input: fermion mass
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // NEW IMPLEMENTATION FIXING FLOATING POINT EXCEPTIONS IN SIMD CODE (#701)
    // Variables xxxDENOM are a hack to avoid division-by-0 FPE while preserving speed (#701 and #727)
    // Variables xxxDENOM are declared as 'volatile' to make sure they are not optimized away on clang! (#724)
    // A few additional variables are declared as 'volatile' to avoid sqrt-of-negative-number FPEs (#736)
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fi = W_ACCESS::kernelAccess( wavefunctions );
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    if( fmass != 0. )
    {
#ifndef MGONGPU_CPPSIMD
      const fptype_sv pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
#else
      volatile fptype_sv p2 = pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3; // volatile fixes #736
      const fptype_sv pp = fpmin( pvec0, fpsqrt( p2 ) );
#endif
      // In C++ ixxxxx, use a single ip/im numbering that is valid both for pp==0 and pp>0, which have two numbering schemes in Fortran ixxxxx:
      // for pp==0, Fortran sqm(0:1) has indexes 0,1 as in C++; but for Fortran pp>0, omega(2) has indexes 1,2 and not 0,1
      // NB: this is only possible in ixxxx, but in oxxxxx two different numbering schemes must be used
      const int ip = ( 1 + nh ) / 2; // NB: same as in Fortran pp==0, differs from Fortran pp>0, which is (3+nh)/2 because omega(2) has indexes 1,2
      const int im = ( 1 - nh ) / 2; // NB: same as in Fortran pp==0, differs from Fortran pp>0, which is (3-nh)/2 because omega(2) has indexes 1,2
#ifndef MGONGPU_CPPSIMD
      if( pp == 0. )
      {
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
        //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
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
        const fptype sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype pp3 = fpmax( pp + pvec3, 0. );
        const cxtype chi[2] = { cxmake( fpsqrt( pp3 * (fptype)0.5 / pp ), 0. ),
                                ( pp3 == 0. ? cxmake( -nh, 0. ) : cxmake( nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        fi[2] = sfomega[0] * chi[im];
        fi[3] = sfomega[0] * chi[ip];
        fi[4] = sfomega[1] * chi[im];
        fi[5] = sfomega[1] * chi[ip];
      }
#else
      // Branch A: pp == 0.
      // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
      fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0 }; // possibility of negative fermion masses (NB: SCALAR!)
      sqm[1] = ( fmass < 0 ? -sqm[0] : sqm[0] );          // AV: removed an abs here (as above)
      const cxtype fiA_2 = ip * sqm[ip];                  // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_3 = im * nsf * sqm[ip];            // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_4 = ip * nsf * sqm[im];            // scalar cxtype: real part initialised from fptype, imag part = 0
      const cxtype fiA_5 = im * sqm[im];                  // scalar cxtype: real part initialised from fptype, imag part = 0
      // Branch B: pp != 0.
      const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                             fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
      fptype_v omega[2] = { fpsqrt( pvec0 + pp ), 0 };
      omega[1] = fmass / omega[0];
      const fptype_v sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const fptype_v pp3 = fpmax( pp + pvec3, 0 );
      volatile fptype_v ppDENOM = fpternary( pp != 0, pp, 1. );    // hack: ppDENOM[ieppV]=1 if pp[ieppV]==0
      volatile fptype_v pp3DENOM = fpternary( pp3 != 0, pp3, 1. ); // hack: pp3DENOM[ieppV]=1 if pp3[ieppV]==0
      volatile fptype_v chi0r2 = pp3 * 0.5 / ppDENOM;              // volatile fixes #736
      const cxtype_v chi[2] = { cxmake( fpsqrt( chi0r2 ), 0 ),     // hack: dummy[ieppV] is not used if pp[ieppV]==0
                                cxternary( ( pp3 == 0. ),
                                           cxmake( -nh, 0 ),
                                           cxmake( (fptype)nh * pvec1, pvec2 ) / fpsqrt( 2. * ppDENOM * pp3DENOM ) ) }; // hack: dummy[ieppV] is not used if pp[ieppV]==0
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
#ifdef MGONGPU_CPPSIMD
      volatile fptype_sv p0p3 = fpmax( pvec0 + pvec3, 0 ); // volatile fixes #736
      volatile fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                             fptype_sv{ 0 },
                                             fpsqrt( p0p3 ) * (fptype)nsf );
      volatile fptype_sv sqp0p3DENOM = fpternary( sqp0p3 != 0, (fptype_sv)sqp0p3, 1. ); // hack: dummy sqp0p3DENOM[ieppV]=1 if sqp0p3[ieppV]==0
      cxtype_sv chi[2] = { cxmake( (fptype_v)sqp0p3, 0. ),
                           cxternary( sqp0p3 == 0,
                                      cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0 ), 0. ),
                                      cxmake( (fptype)nh * pvec1, pvec2 ) / (const fptype_v)sqp0p3DENOM ) }; // hack: dummy[ieppV] is not used if sqp0p3[ieppV]==0
#else
      const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                          fptype_sv{ 0 },
                                          fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
      const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ),
                                 ( sqp0p3 == 0. ? cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0 ), 0. ) : cxmake( (fptype)nh * pvec1, pvec2 ) / sqp0p3 ) };
#endif
      if( nh == 1 )
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
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  ipzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fi = W_ACCESS::kernelAccess( wavefunctions );
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
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  imzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fi = W_ACCESS::kernelAccess( wavefunctions );
    fi[0] = cxmake( pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv chi = cxmake( -(fptype)nhel * fpsqrt( -2. * pvec3 ), 0. );
    fi[3] = cxzero_sv();
    fi[4] = cxzero_sv();
    if( nh == 1 )
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
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  ixzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fi = W_ACCESS::kernelAccess( wavefunctions );
    //fi[0] = cxmake( -pvec0 * nsf, -pvec2 * nsf ); // AV: BUG! not the same as ixxxxx
    //fi[1] = cxmake( -pvec0 * nsf, -pvec1 * nsf ); // AV: BUG! not the same as ixxxxx
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf ); // AV: BUG FIX
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf ); // AV: BUG FIX
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1 / sqp0p3, pvec2 / sqp0p3 );
    if( nh == 1 )
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

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  vxxxxx( const fptype momenta[], // input: momenta
          const fptype vmass,     // input: vector boson mass
          const int nhel,         // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
          const int nsv,          // input: +1 (final) or -1 (initial)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // NEW IMPLEMENTATION FIXING FLOATING POINT EXCEPTIONS IN SIMD CODE (#701)
    // Variables xxxDENOM are a hack to avoid division-by-0 FPE while preserving speed (#701 and #727)
    // Variables xxxDENOM are declared as 'volatile' to make sure they are not optimized away on clang! (#724)
    // A few additional variables are declared as 'volatile' to avoid sqrt-of-negative-number FPEs (#736)
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* vc = W_ACCESS::kernelAccess( wavefunctions );
    const fptype sqh = fpsqrt( 0.5 ); // AV this is > 0!
    const fptype hel = nhel;
    vc[0] = cxmake( pvec0 * (fptype)nsv, pvec3 * (fptype)nsv );
    vc[1] = cxmake( pvec1 * (fptype)nsv, pvec2 * (fptype)nsv );
    if( vmass != 0. )
    {
      const int nsvahl = nsv * std::abs( hel );
      const fptype hel0 = 1. - std::abs( hel );
#ifndef MGONGPU_CPPSIMD
      const fptype_sv pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
      const fptype_sv pp = fpmin( pvec0, fpsqrt( pt2 + ( pvec3 * pvec3 ) ) );
      const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
      if( pp == 0. )
      {
        vc[2] = cxmake( 0., 0. );
        vc[3] = cxmake( -hel * sqh, 0. );
        vc[4] = cxmake( 0., nsvahl * sqh );
        vc[5] = cxmake( hel0, 0. );
      }
      else
      {
        const fptype emp = pvec0 / ( vmass * pp );
        vc[2] = cxmake( hel0 * pp / vmass, 0. );
        vc[5] = cxmake( hel0 * pvec3 * emp + hel * pt / pp * sqh, 0. );
        if( pt != 0. )
        {
          const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
          vc[3] = cxmake( hel0 * pvec1 * emp - pvec1 * pzpt, -nsvahl * pvec2 / pt * sqh );
          vc[4] = cxmake( hel0 * pvec2 * emp - pvec2 * pzpt, nsvahl * pvec1 / pt * sqh );
        }
        else
        {
          vc[3] = cxmake( -hel * sqh, 0. );
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          //vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
          vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
        }
      }
#else
      volatile fptype_sv pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
      volatile fptype_sv p2 = pt2 + ( pvec3 * pvec3 ); // volatile fixes #736
      const fptype_sv pp = fpmin( pvec0, fpsqrt( p2 ) );
      const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
      // Branch A: pp == 0.
      const cxtype vcA_2 = cxmake( 0, 0 );
      const cxtype vcA_3 = cxmake( -hel * sqh, 0 );
      const cxtype vcA_4 = cxmake( 0, nsvahl * sqh );
      const cxtype vcA_5 = cxmake( hel0, 0 );
      // Branch B: pp != 0.
      volatile fptype_v ppDENOM = fpternary( pp != 0, pp, 1. ); // hack: ppDENOM[ieppV]=1 if pp[ieppV]==0
      const fptype_v emp = pvec0 / ( vmass * ppDENOM );         // hack: dummy[ieppV] is not used if pp[ieppV]==0
      const cxtype_v vcB_2 = cxmake( hel0 * pp / vmass, 0 );
      const cxtype_v vcB_5 = cxmake( hel0 * pvec3 * emp + hel * pt / ppDENOM * sqh, 0 ); // hack: dummy[ieppV] is not used if pp[ieppV]==0
      // Branch B1: pp != 0. and pt != 0.
      volatile fptype_v ptDENOM = fpternary( pt != 0, pt, 1. );                                                     // hack: ptDENOM[ieppV]=1 if pt[ieppV]==0
      const fptype_v pzpt = pvec3 / ( ppDENOM * ptDENOM ) * sqh * hel;                                              // hack: dummy[ieppV] is not used if pp[ieppV]==0
      const cxtype_v vcB1_3 = cxmake( hel0 * pvec1 * emp - pvec1 * pzpt, -(fptype)nsvahl * pvec2 / ptDENOM * sqh ); // hack: dummy[ieppV] is not used if pt[ieppV]==0
      const cxtype_v vcB1_4 = cxmake( hel0 * pvec2 * emp - pvec2 * pzpt, (fptype)nsvahl * pvec1 / ptDENOM * sqh );  // hack: dummy[ieppV] is not used if pt[ieppV]==0
      // Branch B2: pp != 0. and pt == 0.
      const cxtype vcB2_3 = cxmake( -hel * sqh, 0. );
      const cxtype_v vcB2_4 = cxmake( 0., (fptype)nsvahl * fpternary( ( pvec3 < 0 ), -sqh, sqh ) ); // AV: removed an abs here
      // Choose between the results from branch A and branch B (and from branch B1 and branch B2)
      const bool_v mask = ( pp == 0. );
      const bool_v maskB = ( pt != 0. );
      vc[2] = cxternary( mask, vcA_2, vcB_2 );
      vc[3] = cxternary( mask, vcA_3, cxternary( maskB, vcB1_3, vcB2_3 ) );
      vc[4] = cxternary( mask, vcA_4, cxternary( maskB, vcB1_4, vcB2_4 ) );
      vc[5] = cxternary( mask, vcA_5, vcB_5 );
#endif
    }
    else
    {
      const fptype_sv& pp = pvec0; // NB: rewrite the following as in Fortran, using pp instead of pvec0
#ifndef MGONGPU_CPPSIMD
      const fptype_sv pt = fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) );
#else
      volatile fptype_sv pt2 = pvec1 * pvec1 + pvec2 * pvec2; // volatile fixes #736
      const fptype_sv pt = fpsqrt( pt2 );
#endif
      vc[2] = cxzero_sv();
      vc[5] = cxmake( hel * pt / pp * sqh, 0. );
#ifndef MGONGPU_CPPSIMD
      if( pt != 0. )
      {
        const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
        vc[3] = cxmake( -pvec1 * pzpt, -nsv * pvec2 / pt * sqh );
        vc[4] = cxmake( -pvec2 * pzpt, nsv * pvec1 / pt * sqh );
      }
      else
      {
        vc[3] = cxmake( -hel * sqh, 0. );
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        //vc[4] = cxmake( 0, nsv * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV why abs here?
        vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
      }
#else
      // Branch A: pt != 0.
      volatile fptype_v ptDENOM = fpternary( pt != 0, pt, 1. );                             // hack: ptDENOM[ieppV]=1 if pt[ieppV]==0
      const fptype_v pzpt = pvec3 / ( pp * ptDENOM ) * sqh * hel;                           // hack: dummy[ieppV] is not used if pt[ieppV]==0
      const cxtype_v vcA_3 = cxmake( -pvec1 * pzpt, -(fptype)nsv * pvec2 / ptDENOM * sqh ); // hack: dummy[ieppV] is not used if pt[ieppV]==0
      const cxtype_v vcA_4 = cxmake( -pvec2 * pzpt, (fptype)nsv * pvec1 / ptDENOM * sqh );  // hack: dummy[ieppV] is not used if pt[ieppV]==0
      // Branch B: pt == 0.
      const cxtype vcB_3 = cxmake( -(fptype)hel * sqh, 0 );
      const cxtype_v vcB_4 = cxmake( 0, (fptype)nsv * fpternary( ( pvec3 < 0 ), -sqh, sqh ) ); // AV: removed an abs here
      // Choose between the results from branch A and branch B
      const bool_v mask = ( pt != 0. );
      vc[3] = cxternary( mask, vcA_3, vcB_3 );
      vc[4] = cxternary( mask, vcA_4, vcB_4 );
#endif
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  sxxxxx( const fptype momenta[], // input: momenta
          //const fptype,                 // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
          //const int,                    // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
          const int nss,          // input: +1 (final) or -1 (initial)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* sc = W_ACCESS::kernelAccess( wavefunctions );
    sc[2] = cxmake( 1 + fptype_sv{ 0 }, 0 );
    sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
    sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  oxxxxx( const fptype momenta[], // input: momenta
          const fptype fmass,     // input: fermion mass
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // NEW IMPLEMENTATION FIXING FLOATING POINT EXCEPTIONS IN SIMD CODE (#701)
    // Variables xxxDENOM are a hack to avoid division-by-0 FPE while preserving speed (#701 and #727)
    // Variables xxxDENOM are declared as 'volatile' to make sure they are not optimized away on clang! (#724)
    // A few additional variables are declared as 'volatile' to avoid sqrt-of-negative-number FPEs (#736)
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fo = W_ACCESS::kernelAccess( wavefunctions );
    fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    if( fmass != 0. )
    {
#ifndef MGONGPU_CPPSIMD
      const fptype_sv pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
      if( pp == 0. )
      {
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
        //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const int ip = -( ( 1 - nh ) / 2 ) * nhel;  // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        const int im = ( 1 + nh ) / 2 * nhel;       // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        fo[2] = cxmake( im * sqm[std::abs( ip )], 0 );
        fo[3] = cxmake( ip * nsf * sqm[std::abs( ip )], 0 );
        fo[4] = cxmake( im * nsf * sqm[std::abs( im )], 0 );
        fo[5] = cxmake( ip * sqm[std::abs( im )], 0 );
      }
      else
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype omega[2] = { fpsqrt( pvec0 + pp ), 0. };
        omega[1] = fmass / omega[0];
        const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sfomeg[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype pp3 = fpmax( pp + pvec3, 0. );
        const cxtype chi[2] = { cxmake( fpsqrt( pp3 * (fptype)0.5 / pp ), 0. ),
                                ( ( pp3 == 0. ) ? cxmake( -nh, 0. )
                                                : cxmake( nh * pvec1, -pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        fo[2] = sfomeg[1] * chi[im];
        fo[3] = sfomeg[1] * chi[ip];
        fo[4] = sfomeg[0] * chi[im];
        fo[5] = sfomeg[0] * chi[ip];
      }
#else
      volatile fptype_sv p2 = pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3; // volatile fixes #736
      const fptype_sv pp = fpmin( pvec0, fpsqrt( p2 ) );
      // Branch A: pp == 0.
      // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
      fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0 }; // possibility of negative fermion masses
      sqm[1] = ( fmass < 0 ? -sqm[0] : sqm[0] );          // AV: removed an abs here (as above)
      const int ipA = -( ( 1 - nh ) / 2 ) * nhel;
      const int imA = ( 1 + nh ) / 2 * nhel;
      const cxtype foA_2 = imA * sqm[std::abs( ipA )];
      const cxtype foA_3 = ipA * nsf * sqm[std::abs( ipA )];
      const cxtype foA_4 = imA * nsf * sqm[std::abs( imA )];
      const cxtype foA_5 = ipA * sqm[std::abs( imA )];
      // Branch B: pp != 0.
      const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                             fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
      fptype_v omega[2] = { fpsqrt( pvec0 + pp ), 0 };
      omega[1] = fmass / omega[0];
      const int ipB = ( 1 + nh ) / 2;
      const int imB = ( 1 - nh ) / 2;
      const fptype_v sfomeg[2] = { sf[0] * omega[ipB], sf[1] * omega[imB] };
      const fptype_v pp3 = fpmax( pp + pvec3, 0. );
      volatile fptype_v ppDENOM = fpternary( pp != 0, pp, 1. );    // hack: ppDENOM[ieppV]=1 if pp[ieppV]==0
      volatile fptype_v pp3DENOM = fpternary( pp3 != 0, pp3, 1. ); // hack: pp3DENOM[ieppV]=1 if pp3[ieppV]==0
      volatile fptype_v chi0r2 = pp3 * 0.5 / ppDENOM;              // volatile fixes #736
      const cxtype_v chi[2] = { cxmake( fpsqrt( chi0r2 ), 0. ),    // hack: dummy[ieppV] is not used if pp[ieppV]==0
                                ( cxternary( ( pp3 == 0. ),
                                             cxmake( -nh, 0. ),
                                             cxmake( (fptype)nh * pvec1, -pvec2 ) / fpsqrt( 2. * ppDENOM * pp3DENOM ) ) ) }; // hack: dummy[ieppV] is not used if pp[ieppV]==0
      const cxtype_v foB_2 = sfomeg[1] * chi[imB];
      const cxtype_v foB_3 = sfomeg[1] * chi[ipB];
      const cxtype_v foB_4 = sfomeg[0] * chi[imB];
      const cxtype_v foB_5 = sfomeg[0] * chi[ipB];
      // Choose between the results from branch A and branch B
      const bool_v mask = ( pp == 0. );
      fo[2] = cxternary( mask, foA_2, foB_2 );
      fo[3] = cxternary( mask, foA_3, foB_3 );
      fo[4] = cxternary( mask, foA_4, foB_4 );
      fo[5] = cxternary( mask, foA_5, foB_5 );
#endif
    }
    else
    {
#ifdef MGONGPU_CPPSIMD
      volatile fptype_sv p0p3 = fpmax( pvec0 + pvec3, 0 ); // volatile fixes #736
      volatile fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                             fptype_sv{ 0 },
                                             fpsqrt( p0p3 ) * (fptype)nsf );
      volatile fptype_v sqp0p3DENOM = fpternary( sqp0p3 != 0, (fptype_sv)sqp0p3, 1. ); // hack: sqp0p3DENOM[ieppV]=1 if sqp0p3[ieppV]==0
      const cxtype_v chi[2] = { cxmake( (fptype_v)sqp0p3, 0. ),
                                cxternary( ( sqp0p3 == 0. ),
                                           cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 ),
                                           cxmake( (fptype)nh * pvec1, -pvec2 ) / (const fptype_sv)sqp0p3DENOM ) }; // hack: dummy[ieppV] is not used if sqp0p3[ieppV]==0
#else
      const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. ) and ( pvec2 == 0. ) and ( pvec3 < 0. ),
                                          0,
                                          fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
      const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ),
                                 ( sqp0p3 == 0. ? cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 ) : cxmake( (fptype)nh * pvec1, -pvec2 ) / sqp0p3 ) };
#endif
      if( nh == 1 )
      {
        fo[2] = chi[0];
        fo[3] = chi[1];
        fo[4] = cxzero_sv();
        fo[5] = cxzero_sv();
      }
      else
      {
        fo[2] = cxzero_sv();
        fo[3] = cxzero_sv();
        fo[4] = chi[1];
        fo[5] = chi[0];
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  opzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fo = W_ACCESS::kernelAccess( wavefunctions );
    fo[0] = cxmake( pvec3 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv csqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
    fo[3] = cxzero_sv();
    fo[4] = cxzero_sv();
    if( nh == 1 )
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
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  omzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fo = W_ACCESS::kernelAccess( wavefunctions );
    fo[0] = cxmake( -pvec3 * (fptype)nsf, pvec3 * (fptype)nsf ); // remember pvec0 == -pvec3
    fo[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv chi1 = cxmake( -nhel, 0. ) * fpsqrt( -2. * pvec3 );
    if( nh == 1 )
    {
      fo[2] = cxzero_sv();
      fo[3] = chi1;
      fo[4] = cxzero_sv();
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[3] = cxzero_sv();
      fo[4] = chi1;
      //fo[5] = chi1; // AV: BUG!
      fo[5] = cxzero_sv(); // AV: BUG FIX
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  template<class M_ACCESS, class W_ACCESS>
  __host__ __device__ void
  oxzxxx( const fptype momenta[], // input: momenta
          //const fptype fmass,   // [skip: ASSUME fermion mass==0]
          const int nhel,         // input: -1 or +1 (helicity of fermion)
          const int nsf,          // input: +1 (particle) or -1 (antiparticle)
          fptype wavefunctions[], // output: wavefunctions
          const int ipar )        // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1 = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2 = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3 = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fo = W_ACCESS::kernelAccess( wavefunctions );
    fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
    if( nh == 1 )
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

  //==========================================================================

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ INLINE void
  VVV1P0_1( const fptype allV2[],
            const fptype allV3[],
            const fptype allCOUP[],
            const double Ccoeff,
            const fptype M1,
            const fptype W1,
            fptype allV1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  template<class W_ACCESS, class A_ACCESS, class C_ACCESS>
  __device__ INLINE void
  FFV1_0( const fptype allF1[],
          const fptype allF2[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          fptype allvertexes[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ INLINE void
  FFV1_1( const fptype allF2[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          const fptype M1,
          const fptype W1,
          fptype allF1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ INLINE void
  FFV1_2( const fptype allF1[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          const fptype M2,
          const fptype W2,
          fptype allF2[] ) ALWAYS_INLINE;

  //==========================================================================

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ void
  VVV1P0_1( const fptype allV2[],
            const fptype allV3[],
            const fptype allCOUP[],
            const double Ccoeff,
            const fptype M1,
            const fptype W1,
            fptype allV1[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv* V2 = W_ACCESS::kernelAccessConst( allV2 );
    const cxtype_sv* V3 = W_ACCESS::kernelAccessConst( allV3 );
    const cxtype_sv COUP = C_ACCESS::kernelAccessConst( allCOUP );
    cxtype_sv* V1 = W_ACCESS::kernelAccess( allV1 );
    const cxtype cI = cxmake( 0., 1. );
    const fptype_sv P2[4] = { +cxreal( V2[0] ), +cxreal( V2[1] ), +cximag( V2[1] ), +cximag( V2[0] ) };
    const fptype_sv P3[4] = { +cxreal( V3[0] ), +cxreal( V3[1] ), +cximag( V3[1] ), +cximag( V3[0] ) };
    V1[0] = +V2[0] + V3[0];
    V1[1] = +V2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP1 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP2 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP3 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP4 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( ( P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( TMP4 * ( -cI * P2[0] + cI * P3[0] ) + ( V2[2] * ( -cI * TMP0 + cI * TMP1 ) + V3[2] * ( +cI * TMP2 - cI * TMP3 ) ) );
    V1[3] = denom * ( TMP4 * ( -cI * P2[1] + cI * P3[1] ) + ( V2[3] * ( -cI * TMP0 + cI * TMP1 ) + V3[3] * ( +cI * TMP2 - cI * TMP3 ) ) );
    V1[4] = denom * ( TMP4 * ( -cI * P2[2] + cI * P3[2] ) + ( V2[4] * ( -cI * TMP0 + cI * TMP1 ) + V3[4] * ( +cI * TMP2 - cI * TMP3 ) ) );
    V1[5] = denom * ( TMP4 * ( -cI * P2[3] + cI * P3[3] ) + ( V2[5] * ( -cI * TMP0 + cI * TMP1 ) + V3[5] * ( +cI * TMP2 - cI * TMP3 ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  template<class W_ACCESS, class A_ACCESS, class C_ACCESS>
  __device__ void
  FFV1_0( const fptype allF1[],
          const fptype allF2[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          fptype allvertexes[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv* F1 = W_ACCESS::kernelAccessConst( allF1 );
    const cxtype_sv* F2 = W_ACCESS::kernelAccessConst( allF2 );
    const cxtype_sv* V3 = W_ACCESS::kernelAccessConst( allV3 );
    const cxtype_sv COUP = C_ACCESS::kernelAccessConst( allCOUP );
    cxtype_sv* vertex = A_ACCESS::kernelAccess( allvertexes );
    const cxtype cI = cxmake( 0., 1. );
    const cxtype_sv TMP5 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) + ( F1[3] * ( F2[4] * ( V3[3] - cI * V3[4] ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * V3[4] ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * V3[4] ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    ( *vertex ) = COUP * -cI * TMP5;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ void
  FFV1_1( const fptype allF2[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          const fptype M1,
          const fptype W1,
          fptype allF1[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv* F2 = W_ACCESS::kernelAccessConst( allF2 );
    const cxtype_sv* V3 = W_ACCESS::kernelAccessConst( allV3 );
    const cxtype_sv COUP = C_ACCESS::kernelAccessConst( allCOUP );
    cxtype_sv* F1 = W_ACCESS::kernelAccess( allF1 );
    const cxtype cI = cxmake( 0., 1. );
    F1[0] = +F2[0] + V3[0];
    F1[1] = +F2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( F1[0] ), -cxreal( F1[1] ), -cximag( F1[1] ), -cximag( F1[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( ( P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    F1[2] = denom * cI * ( F2[2] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3] - cI * V3[4] ) + ( P1[2] * ( +cI * V3[3] + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[3] + cI * V3[4] ) + ( P1[1] * ( -one ) * ( V3[2] + V3[5] ) + ( P1[2] * ( -one ) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( V3[3] + cI * V3[4] ) ) ) ) + M1 * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) ) );
    F1[3] = denom * ( -cI ) * ( F2[2] * ( P1[0] * ( -V3[3] + cI * V3[4] ) + ( P1[1] * ( V3[2] - V3[5] ) + ( P1[2] * ( -cI * V3[2] + cI * V3[5] ) + P1[3] * ( V3[3] - cI * V3[4] ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * ( -one ) * ( V3[3] + cI * V3[4] ) + ( P1[2] * ( +cI * V3[3] - V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + M1 * ( F2[4] * ( -V3[3] + cI * V3[4] ) + F2[5] * ( -V3[2] + V3[5] ) ) ) );
    F1[4] = denom * ( -cI ) * ( F2[4] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * ( -V3[3] + cI * V3[4] ) + ( P1[2] * ( -one ) * ( +cI * V3[3] + V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + ( F2[5] * ( P1[0] * ( V3[3] + cI * V3[4] ) + ( P1[1] * ( -V3[2] + V3[5] ) + ( P1[2] * ( -cI * V3[2] + cI * V3[5] ) - P1[3] * ( V3[3] + cI * V3[4] ) ) ) ) + M1 * ( F2[2] * ( -V3[2] + V3[5] ) + F2[3] * ( V3[3] + cI * V3[4] ) ) ) );
    F1[5] = denom * cI * ( F2[4] * ( P1[0] * ( -V3[3] + cI * V3[4] ) + ( P1[1] * ( V3[2] + V3[5] ) + ( P1[2] * ( -one ) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( -V3[3] + cI * V3[4] ) ) ) ) + ( F2[5] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3] + cI * V3[4] ) + ( P1[2] * ( -cI * V3[3] + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + M1 * ( F2[2] * ( -V3[3] + cI * V3[4] ) + F2[3] * ( V3[2] + V3[5] ) ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ void
  FFV1_2( const fptype allF1[],
          const fptype allV3[],
          const fptype allCOUP[],
          const double Ccoeff,
          const fptype M2,
          const fptype W2,
          fptype allF2[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv* F1 = W_ACCESS::kernelAccessConst( allF1 );
    const cxtype_sv* V3 = W_ACCESS::kernelAccessConst( allV3 );
    const cxtype_sv COUP = C_ACCESS::kernelAccessConst( allCOUP );
    cxtype_sv* F2 = W_ACCESS::kernelAccess( allF2 );
    const cxtype cI = cxmake( 0., 1. );
    F2[0] = +F1[0] + V3[0];
    F2[1] = +F1[1] + V3[1];
    const fptype_sv P2[4] = { -cxreal( F2[0] ), -cxreal( F2[1] ), -cximag( F2[1] ), -cximag( F2[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( ( P2[0] * P2[0] ) - ( P2[1] * P2[1] ) - ( P2[2] * P2[2] ) - ( P2[3] * P2[3] ) - M2 * ( M2 - cI * W2 ) );
    F2[2] = denom * cI * ( F1[2] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * ( -one ) * ( V3[3] + cI * V3[4] ) + ( P2[2] * ( +cI * V3[3] - V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + ( F1[3] * ( P2[0] * ( V3[3] - cI * V3[4] ) + ( P2[1] * ( -V3[2] + V3[5] ) + ( P2[2] * ( +cI * V3[2] - cI * V3[5] ) + P2[3] * ( -V3[3] + cI * V3[4] ) ) ) ) + M2 * ( F1[4] * ( V3[2] - V3[5] ) + F1[5] * ( -V3[3] + cI * V3[4] ) ) ) );
    F2[3] = denom * ( -cI ) * ( F1[2] * ( P2[0] * ( -one ) * ( V3[3] + cI * V3[4] ) + ( P2[1] * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) - P2[3] * ( V3[3] + cI * V3[4] ) ) ) ) + ( F1[3] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3] - cI * V3[4] ) + ( P2[2] * ( +cI * V3[3] + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + M2 * ( F1[4] * ( V3[3] + cI * V3[4] ) - F1[5] * ( V3[2] + V3[5] ) ) ) );
    F2[4] = denom * ( -cI ) * ( F1[4] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3] + cI * V3[4] ) + ( P2[2] * ( -cI * V3[3] + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[3] - cI * V3[4] ) + ( P2[1] * ( -one ) * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) + P2[3] * ( V3[3] - cI * V3[4] ) ) ) ) + M2 * ( F1[2] * ( -one ) * ( V3[2] + V3[5] ) + F1[3] * ( -V3[3] + cI * V3[4] ) ) ) );
    F2[5] = denom * cI * ( F1[4] * ( P2[0] * ( -one ) * ( V3[3] + cI * V3[4] ) + ( P2[1] * ( V3[2] - V3[5] ) + ( P2[2] * ( +cI * V3[2] - cI * V3[5] ) + P2[3] * ( V3[3] + cI * V3[4] ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * ( -V3[3] + cI * V3[4] ) + ( P2[2] * ( -one ) * ( +cI * V3[3] + V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + M2 * ( F1[2] * ( V3[3] + cI * V3[4] ) + F1[3] * ( V3[2] - V3[5] ) ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace

#endif // HelAmps_sm_H
