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
// MadGraph5_aMC@NLO v. 3.5.0_lo_vect, 2023-06-09
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
    // USE _SV SUFFIXES EVERYWHERE AND ITERATE EXPLICITLY OVER ALL ELEMENTS OF SIMD VECTORS IF NEEDED
    const fptype_sv& pvec0_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fi_sv = W_ACCESS::kernelAccess( wavefunctions );
    fi_sv[0] = cxmake( -pvec0_sv * (fptype)nsf, -pvec3_sv * (fptype)nsf );
    fi_sv[1] = cxmake( -pvec1_sv * (fptype)nsf, -pvec2_sv * (fptype)nsf );
    const int nh = nhel * nsf;
    if( fmass != 0. )
    {
      const fptype_sv pp_sv = fpmin( pvec0_sv, fpsqrt( pvec1_sv * pvec1_sv + pvec2_sv * pvec2_sv + pvec3_sv * pvec3_sv ) );
      const fptype_sv pp3_sv = fpmax( pp_sv + pvec3_sv, 0. );
      // In C++ ixxxxx, use a single ip/im numbering that is valid both for pp==0 and pp>0, which have two numbering schemes in Fortran ixxxxx:
      // for pp==0, Fortran sqm(0:1) has indexes 0,1 as in C++; but for Fortran pp>0, omega(2) has indexes 1,2 and not 0,1
      // NB: this is only possible in ixxxx, but in oxxxxx two different numbering schemes must be used
      const int ip = ( 1 + nh ) / 2; // NB: same as in Fortran pp==0, differs from Fortran pp>0, which is (3+nh)/2 because omega(2) has indexes 1,2
      const int im = ( 1 - nh ) / 2; // NB: same as in Fortran pp==0, differs from Fortran pp>0, which is (3-nh)/2 because omega(2) has indexes 1,2
      if( maskand( pp_sv == 0. ) )
      {
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs around sqm[0] here
        fi_sv[2] = cxtype_sv( ip * sqm[ip] );       // IIII=0000
        fi_sv[3] = cxtype_sv( im * nsf * sqm[ip] ); // IIII=0000
        fi_sv[4] = cxtype_sv( ip * nsf * sqm[im] ); // IIII=0000
        fi_sv[5] = cxtype_sv( im * sqm[im] );       // IIII=0000
      }
      else if( maskand( pp_sv != 0. ) and maskand( pp3_sv == 0. ) )
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        const fptype_sv sfomega_sv[2] = { sf[0] * omega_sv[ip], sf[1] * omega_sv[im] };
        const cxtype_sv chi_sv[2] = { cxmake( fpsqrt( pp3_sv * (fptype)0.5 / pp_sv ), 0. ),
                                      cxtype_sv( -nh ) }; // IIII=0000
        fi_sv[2] = sfomega_sv[0] * chi_sv[im];
        fi_sv[3] = sfomega_sv[0] * chi_sv[ip];
        fi_sv[4] = sfomega_sv[1] * chi_sv[im];
        fi_sv[5] = sfomega_sv[1] * chi_sv[ip];
      }
      else if( maskand( pp_sv != 0. ) and maskand( pp3_sv != 0. ) )
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        const fptype_sv sfomega_sv[2] = { sf[0] * omega_sv[ip], sf[1] * omega_sv[im] };
        const cxtype_sv chi_sv[2] = { cxmake( fpsqrt( pp3_sv * (fptype)0.5 / pp_sv ), 0. ),
                                      cxmake( (fptype)nh * pvec1_sv, pvec2_sv ) / fpsqrt( 2. * pp_sv * pp3_sv ) };
        fi_sv[2] = sfomega_sv[0] * chi_sv[im];
        fi_sv[3] = sfomega_sv[0] * chi_sv[ip];
        fi_sv[4] = sfomega_sv[1] * chi_sv[im];
        fi_sv[5] = sfomega_sv[1] * chi_sv[ip];
      }
      else
      {
#ifdef MGONGPU_CPPSIMD
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses (NB: use std::abs!)
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        const fptype_sv sfomega_sv[2] = { sf[0] * omega_sv[ip], sf[1] * omega_sv[im] };
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
        {
          const fptype& pp = pp_sv[ieppV];
          if( pp == 0. )
          {
            fi_sv[2][ieppV] = cxmake( ip * sqm[ip], 0 );
            fi_sv[3][ieppV] = cxmake( im * nsf * sqm[ip], 0 );
            fi_sv[4][ieppV] = cxmake( ip * nsf * sqm[im], 0 );
            fi_sv[5][ieppV] = cxmake( im * sqm[im], 0 );
          }
          else
          {
            const fptype& pp3 = pp3_sv[ieppV];
            const fptype& pvec1 = pvec1_sv[ieppV];
            const fptype& pvec2 = pvec2_sv[ieppV];
            const cxtype chi[2] = { cxmake( fpsqrt( pp3 * (fptype)0.5 / pp ), 0. ),
                                    ( pp3 == 0. ? cxmake( -nh, 0. ) : cxmake( nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
            fi_sv[2][ieppV] = sfomega_sv[0][ieppV] * chi[im];
            fi_sv[3][ieppV] = sfomega_sv[0][ieppV] * chi[ip];
            fi_sv[4][ieppV] = sfomega_sv[1][ieppV] * chi[im];
            fi_sv[5][ieppV] = sfomega_sv[1][ieppV] * chi[ip];
          }
        }
#else
        printf( "INTERNAL ERROR in ixxxxx: no path to this statement on GPUs or scalar C++!\n" );
        assert( false );
#endif
      }
    }
    else
    {
      const fptype_sv sqp0p3_sv = fpternary( ( pvec1_sv == 0. and pvec2_sv == 0. and pvec3_sv < 0. ),
                                             fptype_sv{ 0 },
                                             fpsqrt( fpmax( pvec0_sv + pvec3_sv, 0. ) ) * (fptype)nsf );
#ifdef MGONGPU_CPPSIMD
      cxtype_sv chi_sv[2] = { cxmake( sqp0p3_sv, 0. ), cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0_sv ), 0. ) };
      std::cout << "Entering loop" << std::endl;
      // This loop triggers a FPE in no-debug builds (most likely because it is auto-vectorized)
#pragma GCC push_options
#pragma GCC optimize("O0")
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        if( sqp0p3_sv[ieppV] != 0. )
        {
          chi_sv[1][ieppV] = cxmake( (fptype)nh * pvec1_sv[ieppV], pvec2_sv[ieppV] ) / sqp0p3_sv[ieppV];
        }
      }
#pragma GCC pop_options
      std::cout << "Completed loop" << std::endl;
#else
      const cxtype_sv chi_sv[2] = { cxmake( sqp0p3_sv, 0. ),
                                    ( sqp0p3_sv == 0. ? cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0_sv ), 0. ) : cxmake( (fptype)nh * pvec1_sv, pvec2_sv ) / sqp0p3_sv ) };
#endif
      if( nh == 1 )
      {
        fi_sv[2] = cxzero_sv();
        fi_sv[3] = cxzero_sv();
        fi_sv[4] = chi_sv[0];
        fi_sv[5] = chi_sv[1];
      }
      else
      {
        fi_sv[2] = chi_sv[1];
        fi_sv[3] = chi_sv[0];
        fi_sv[4] = cxzero_sv();
        fi_sv[5] = cxzero_sv();
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
    // USE _SV SUFFIXES EVERYWHERE AND ITERATE EXPLICITLY OVER ALL ELEMENTS OF SIMD VECTORS IF NEEDED
    const fptype_sv& pvec0_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* vc_sv = W_ACCESS::kernelAccess( wavefunctions );
    const fptype sqh = fpsqrt( 0.5 ); // AV this is > 0!
    const fptype hel = nhel;
    vc_sv[0] = cxmake( pvec0_sv * (fptype)nsv, pvec3_sv * (fptype)nsv );
    vc_sv[1] = cxmake( pvec1_sv * (fptype)nsv, pvec2_sv * (fptype)nsv );
    if( vmass != 0. )
    {
      const int nsvahl = nsv * std::abs( hel );
      const fptype_sv pt2_sv = ( pvec1_sv * pvec1_sv ) + ( pvec2_sv * pvec2_sv );
      const fptype_sv pp_sv = fpmin( pvec0_sv, fpsqrt( pt2_sv + ( pvec3_sv * pvec3_sv ) ) );
      const fptype_sv pt_sv = fpmin( pp_sv, fpsqrt( pt2_sv ) );
      const fptype hel0 = 1. - std::abs( hel );
      if( maskand( pp_sv == 0. ) )
      {
        vc_sv[2] = cxzero_sv();
        vc_sv[3] = cxtype_sv( -hel * sqh );                                // IIII=0000
        vc_sv[4] = cxtype_sv( fptype_sv{}, fptype_sv{} + (nsvahl * sqh) ); // RRRR=0000
        vc_sv[5] = cxtype_sv( hel0 );                                      // IIII=0000
      }
      else if( maskand( pp_sv != 0. ) && maskand( pt_sv == 0. ) )
      {
        const fptype_sv emp_sv = pvec0_sv / ( vmass * pp_sv );
        vc_sv[2] = cxtype_sv( hel0 * pp_sv / vmass );                                                // IIII=0000
        vc_sv[5] = cxtype_sv( hel0 * pvec3_sv * emp_sv + hel * pt_sv / pp_sv * sqh );                // IIII=0000
        vc_sv[3] = cxtype_sv( -hel * sqh );                                                          // IIII=0000
        vc_sv[4] = cxtype_sv( fptype_sv{}, (fptype)nsvahl * fpternary( pvec3_sv < 0., -sqh, sqh ) ); // AV: removed an abs here
      }
      else if( maskand( pp_sv != 0. ) && maskand( pt_sv != 0. ) )
      {
        const fptype_sv emp_sv = pvec0_sv / ( vmass * pp_sv );
        vc_sv[2] = cxtype_sv( hel0 * pp_sv / vmass );                                 // IIII=0000
        vc_sv[5] = cxtype_sv( hel0 * pvec3_sv * emp_sv + hel * pt_sv / pp_sv * sqh ); // IIII=0000
        const fptype_sv pzpt_sv = pvec3_sv / ( pp_sv * pt_sv ) * sqh * hel;
        vc_sv[3] = cxmake( hel0 * pvec1_sv * emp_sv - pvec1_sv * pzpt_sv, -(fptype)nsvahl * pvec2_sv / pt_sv * sqh );
        vc_sv[4] = cxmake( hel0 * pvec2_sv * emp_sv - pvec2_sv * pzpt_sv, (fptype)nsvahl * pvec1_sv / pt_sv * sqh );
      }
      else
      {
#ifdef MGONGPU_CPPSIMD
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
        {
          const fptype& pp = pp_sv[ieppV];
          const fptype& pvec0 = pvec0_sv[ieppV];
          const fptype& pvec1 = pvec1_sv[ieppV];
          const fptype& pvec2 = pvec2_sv[ieppV];
          const fptype& pvec3 = pvec3_sv[ieppV];
          if( pp == 0. )
          {
            vc_sv[2][ieppV] = cxtype( 0., 0. );
            vc_sv[3][ieppV] = cxtype( -hel * sqh, 0. );
            vc_sv[4][ieppV] = cxtype( 0., nsvahl * sqh );
            vc_sv[5][ieppV] = cxtype( hel0, 0. );
          }
          else
          {
            const fptype& pt = pt_sv[ieppV];
            const fptype emp = pvec0 / ( vmass * pp );
            vc_sv[2][ieppV] = cxtype( hel0 * pp / vmass, 0. );
            vc_sv[5][ieppV] = cxtype( hel0 * pvec3 * emp + hel * pt / pp * sqh );
            if ( pt == 0. )
            {
              vc_sv[3][ieppV] = cxtype( -hel * sqh, 0. );
              vc_sv[4][ieppV] = cxtype( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
            }
            else
            {
              const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
              vc_sv[3][ieppV] = cxtype( hel0 * pvec1 * emp - pvec1 * pzpt, -nsvahl * pvec2 / pt * sqh );
              vc_sv[4][ieppV] = cxtype( hel0 * pvec2 * emp - pvec2 * pzpt, nsvahl * pvec1 / pt * sqh );
            }
          }
        }        
#else
        printf( "INTERNAL ERROR in ixxxxx: no path to this statement on GPUs or scalar C++!\n" );
        assert( false );
#endif
      }
    }
    else
    {
      const fptype_sv& pp_sv = pvec0_sv; // NB: rewrite the following as in Fortran, using pp instead of pvec0
      const fptype_sv pt_sv = fpsqrt( ( pvec1_sv * pvec1_sv ) + ( pvec2_sv * pvec2_sv ) );
      vc_sv[2] = cxzero_sv();
      vc_sv[5] = cxmake( hel * pt_sv / pp_sv * sqh, 0. );
      if( maskand( pt_sv == 0. ) )
      {
        vc_sv[3] = cxtype_sv( -hel * sqh );                                                       // IIII=0000
        vc_sv[4] = cxtype_sv( fptype_sv{}, (fptype)nsv * fpternary( pvec3_sv < 0., -sqh, sqh ) ); // AV: removed an abs here
      }
      else if( maskand( pt_sv != 0. ) )
      {
        const fptype_sv pzpt_sv = pvec3_sv / ( pp_sv * pt_sv ) * sqh * hel;
        vc_sv[3] = cxmake( -pvec1_sv * pzpt_sv, -(fptype)nsv * pvec2_sv / pt_sv * sqh );
        vc_sv[4] = cxmake( -pvec2_sv * pzpt_sv, (fptype)nsv * pvec1_sv / pt_sv * sqh );
      }        
      else
      {
#ifdef MGONGPU_CPPSIMD
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
        {
          const fptype& pp = pp_sv[ieppV];
          const fptype& pt = pt_sv[ieppV];
          const fptype& pvec1 = pvec1_sv[ieppV];
          const fptype& pvec2 = pvec2_sv[ieppV];
          const fptype& pvec3 = pvec3_sv[ieppV];
          if( pt == 0. )
          {
            vc_sv[3][ieppV] = cxtype( -hel * sqh );
            vc_sv[4][ieppV] = cxtype( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
          }
          else
          {
            const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
            vc_sv[3][ieppV] = cxtype( -pvec1 * pzpt, -nsv * pvec2 / pt * sqh );
            vc_sv[4][ieppV] = cxtype( -pvec2 * pzpt, nsv * pvec1 / pt * sqh );
          }
        }
#else
        printf( "INTERNAL ERROR in vxxxxx: no path to this statement on GPUs or scalar C++!\n" );
        assert( false );
#endif
      }
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
    // USE _SV SUFFIXES EVERYWHERE AND ITERATE EXPLICITLY OVER ALL ELEMENTS OF SIMD VECTORS IF NEEDED
    const fptype_sv& pvec0_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 0, ipar );
    const fptype_sv& pvec1_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 1, ipar );
    const fptype_sv& pvec2_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 2, ipar );
    const fptype_sv& pvec3_sv = M_ACCESS::kernelAccessIp4IparConst( momenta, 3, ipar );
    cxtype_sv* fo_sv = W_ACCESS::kernelAccess( wavefunctions );
    fo_sv[0] = cxmake( pvec0_sv * (fptype)nsf, pvec3_sv * (fptype)nsf );
    fo_sv[1] = cxmake( pvec1_sv * (fptype)nsf, pvec2_sv * (fptype)nsf );
    const int nh = nhel * nsf;
    if( fmass != 0. )
    {
      const fptype_sv pp_sv = fpmin( pvec0_sv, fpsqrt( pvec1_sv * pvec1_sv + pvec2_sv * pvec2_sv + pvec3_sv * pvec3_sv ) );
      const fptype_sv pp3_sv = fpmax( pp_sv + pvec3_sv, 0. );
      if( maskand( pp_sv == 0. ) )
      {
        const int ip = -( ( 1 - nh ) / 2 ) * nhel;  // NB: same as in Fortran! Fortran sqm(0:1) also has indexes 0,1 as in C++
        const int im = ( 1 + nh ) / 2 * nhel;       // NB: same as in Fortran! Fortran sqm(0:1) also has indexes 0,1 as in C++
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs around sqm[0] here
        fo_sv[2] = cxtype_sv( im * sqm[std::abs( ip )] );       // IIII=0000
        fo_sv[3] = cxtype_sv( ip * nsf * sqm[std::abs( ip )] ); // IIII=0000
        fo_sv[4] = cxtype_sv( im * nsf * sqm[std::abs( im )] ); // IIII=0000
        fo_sv[5] = cxtype_sv( ip * sqm[std::abs( im )] );       // IIII=0000
      }
      else if( maskand( pp_sv != 0. ) and maskand( pp3_sv == 0. ) )
      {
        const int ip = ( 1 + nh ) / 2; // NB: differs from Fortran! Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: differs from Fortran! Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        const fptype_sv sfomega_sv[2] = { sf[0] * omega_sv[ip], sf[1] * omega_sv[im] };
        const cxtype_sv chi_sv[2] = { cxmake( fpsqrt( pp3_sv * (fptype)0.5 / pp_sv ), 0. ),
                                      cxtype_sv( -nh ) }; // IIII=0000
        fo_sv[2] = sfomega_sv[1] * chi_sv[im];
        fo_sv[3] = sfomega_sv[1] * chi_sv[ip];
        fo_sv[4] = sfomega_sv[0] * chi_sv[im];
        fo_sv[5] = sfomega_sv[0] * chi_sv[ip];
      }
      else if( maskand( pp_sv != 0. ) and maskand( pp3_sv != 0. ) )
      {
        const int ip = ( 1 + nh ) / 2; // NB: differs from Fortran! Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: differs from Fortran! Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        const fptype_sv sfomega_sv[2] = { sf[0] * omega_sv[ip], sf[1] * omega_sv[im] };
        const cxtype_sv chi_sv[2] = { cxmake( fpsqrt( pp3_sv * (fptype)0.5 / pp_sv ), 0. ),
                                      cxmake( (fptype)nh * pvec1_sv, -pvec2_sv ) / fpsqrt( 2. * pp_sv * pp3_sv ) };
        fo_sv[2] = sfomega_sv[1] * chi_sv[im];
        fo_sv[3] = sfomega_sv[1] * chi_sv[ip];
        fo_sv[4] = sfomega_sv[0] * chi_sv[im];
        fo_sv[5] = sfomega_sv[0] * chi_sv[ip];
      }
      else
      {
#ifdef MGONGPU_CPPSIMD
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses (NB: use std::abs!)
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_sv omega_sv[2] = { fpsqrt( pvec0_sv + pp_sv ), 0. };
        omega_sv[1] = fmass / omega_sv[0];
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
        {
          const fptype& pp = pp_sv[ieppV];
          if( pp == 0. )
          {
            const int ip = -( ( 1 - nh ) / 2 ) * nhel;  // NB: same as in Fortran! Fortran sqm(0:1) also has indexes 0,1 as in C++
            const int im = ( 1 + nh ) / 2 * nhel;       // NB: same as in Fortran! Fortran sqm(0:1) also has indexes 0,1 as in C++
            fo_sv[2][ieppV] = cxtype( im * sqm[std::abs( ip )], 0 );
            fo_sv[3][ieppV] = cxtype( ip * nsf * sqm[std::abs( ip )], 0 );
            fo_sv[4][ieppV] = cxtype( im * nsf * sqm[std::abs( im )], 0 );
            fo_sv[5][ieppV] = cxtype( ip * sqm[std::abs( im )], 0 );
          }
          else
          {
            const int ip = ( 1 + nh ) / 2; // NB: differs from Fortran! Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
            const int im = ( 1 - nh ) / 2; // NB: differs from Fortran! Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
            const fptype sfomega[2] = { sf[0] * omega_sv[ip][ieppV], sf[1] * omega_sv[im][ieppV] };
            const fptype& pp3 = pp3_sv[ieppV];
            const fptype& pvec1 = pvec1_sv[ieppV];
            const fptype& pvec2 = pvec2_sv[ieppV];
            const cxtype chi[2] = { cxmake( fpsqrt( pp3 * (fptype)0.5 / pp ), 0. ),
                                    ( pp3 == 0. ? cxmake( -nh, 0. ) : cxmake( nh * pvec1, -pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
            fo_sv[2][ieppV] = sfomega[1] * chi[im];
            fo_sv[3][ieppV] = sfomega[1] * chi[ip];
            fo_sv[4][ieppV] = sfomega[0] * chi[im];
            fo_sv[5][ieppV] = sfomega[0] * chi[ip];
          }
        }
#else
        printf( "INTERNAL ERROR in oxxxxx: no path to this statement on GPUs or scalar C++!\n" );
        assert( false );
#endif
      }
    }
    else
    {
      const fptype_sv sqp0p3_sv = fpternary( ( pvec1_sv == 0. ) and ( pvec2_sv == 0. ) and ( pvec3_sv < 0. ),
                                             fptype_sv{ 0 },
                                             fpsqrt( fpmax( pvec0_sv + pvec3_sv, 0. ) ) * (fptype)nsf );
#ifdef MGONGPU_CPPSIMD
      cxtype_sv chi_sv[2] = { cxmake( sqp0p3_sv, 0. ), cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0_sv ), 0. ) };
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        if( sqp0p3_sv[ieppV] != 0. )
          chi_sv[1][ieppV] = cxmake( (fptype)nh * pvec1_sv[ieppV], -pvec2_sv[ieppV] ) / sqp0p3_sv[ieppV];
      }
#else
      const cxtype_sv chi_sv[2] = { cxmake( sqp0p3_sv, 0. ),
                                    ( sqp0p3_sv == 0. ? cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0_sv ), 0. ) : cxmake( (fptype)nh * pvec1_sv, -pvec2_sv ) / sqp0p3_sv ) };
#endif
      if( nh == 1 )
      {
        fo_sv[2] = chi_sv[0];
        fo_sv[3] = chi_sv[1];
        fo_sv[4] = cxzero_sv();
        fo_sv[5] = cxzero_sv();
      }
      else
      {
        fo_sv[2] = cxzero_sv();
        fo_sv[3] = cxzero_sv();
        fo_sv[4] = chi_sv[1];
        fo_sv[5] = chi_sv[0];
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
          fptype allvertexes[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  template<class W_ACCESS, class C_ACCESS>
  __device__ INLINE void
  FFV1_1( const fptype allF2[],
          const fptype allV3[],
          const fptype allCOUP[],
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
