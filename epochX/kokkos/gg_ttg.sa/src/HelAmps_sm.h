//==========================================================================
// This file has been automatically generated for KOKKOS standalone
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include "mgOnGpuTypes.h"
#include "mgOnGpuConfig.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm {



  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void ixxxxx(const T& pvec,
                            const fptype fmass, 
                            const int nhel,
                            const int nsf,
                            cxtype_sv* fi)
{
  
  
  fi[0] = cxmake(-pvec(0) * nsf, -pvec(3) * nsf);
  fi[1] = cxmake(-pvec(1) * nsf, -pvec(2) * nsf);
  const int nh = nhel * nsf;
  if ( fmass != 0. )
  {
    const fptype_sv pp = fpmin(pvec(0), fpsqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) ) );
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
      fptype omega[2] = { fpsqrt( pvec(0) + pp ), 0. };
      omega[1] = fmass / omega[0];
      const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
      const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
      const fptype sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const fptype pp3 = fpmax( pp + pvec(3), 0. );
      const cxtype chi[2] = { cxmake( fpsqrt ( pp3 * (fptype)0.5 / pp ), 0. ),
                              ( pp3 == 0. ?
                                cxmake( -nh, 0. ) :
                                cxmake( nh * pvec(1), pvec(2) ) / fpsqrt( 2. * pp * pp3 ) ) };
      fi[2] = sfomega[0] * chi[im];
      fi[3] = sfomega[0] * chi[ip];
      fi[4] = sfomega[1] * chi[im];
      fi[5] = sfomega[1] * chi[ip];
    }
    
  }
  else
  { 
    const fptype_sv sqp0p3 = fpternary( ( pvec(1) == 0. and pvec(2) == 0. and pvec(3) < 0. ),
                                        fptype_sv{0}, fpsqrt(fpmax(pvec(0) + pvec(3), 0. ) ) * (fptype)nsf );

    const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ), cxternary( ( sqp0p3 == 0. ),
                                                               cxmake( -(fptype)nhel * fpsqrt( 2. * pvec(0) ), 0. ),
                                                               cxmake( (fptype)nh * pvec(1), pvec(2)) / sqp0p3 ) };
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void ipzxxx(const T& pvec,
                            const int& nhel,
                            const int& nsf,
                            cxtype_sv fi[])
{
  const fptype& pvec3 = pvec(3);
  
  fi[0] = cxmake(-pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void imzxxx(const T& pvec, 
                            const int nhel,
                            const int nsf,
                            cxtype_sv fi[])
{

  const fptype pvec3 = pvec(3);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void ixzxxx(const T& pvec,
                            const int& nhel,
                            const int& nsf,
                            cxtype_sv fi[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);

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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void vxxxxx(const T& pvec,
                            const fptype vmass,
                            const int nhel,
                            const int nsv,
                            cxtype_sv vc[]) 
{
  const fptype sqh = fpsqrt( 0.5 );
  const fptype hel = nhel;
  const int nsvahl = nsv * abs(hel);
  vc[0] = cxmake( pvec(0) * (fptype)nsv, pvec(3) * (fptype)nsv );
  vc[1] = cxmake( pvec(1) * (fptype)nsv, pvec(2) * (fptype)nsv );
  if ( vmass != 0. )
  {
    const int nsvahl = nsv * std::abs( hel );
    const fptype_sv pt2 = ( pvec(1) * pvec(1) ) + ( pvec(2) * pvec(2) );
    const fptype_sv pp = fpmin( pvec(0), fpsqrt( pt2 + ( pvec(3) * pvec(3) ) ) );
    const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
    const fptype hel0 = 1. - std::abs( hel );
    if ( pp == 0. )
    {
      vc[2] = cxmake( 0., 0. );
      vc[3] = cxmake( -hel * sqh, 0. );
      vc[4] = cxmake( 0., nsvahl * sqh );
      vc[5] = cxmake( hel0, 0. );
    }
    else
    {
      const fptype emp = pvec(0) / ( vmass * pp );
      vc[2] = cxmake( hel0 * pp / vmass, 0. );
      vc[5] = cxmake( hel0 * pvec(3) * emp + hel * pt / pp * sqh, 0. );
      if ( pt != 0. )
      {
        const fptype pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
        vc[3] = cxmake( hel0 * pvec(1) * emp - pvec(1) * pzpt, -nsvahl * pvec(2) / pt * sqh );
        vc[4] = cxmake( hel0 * pvec(2) * emp - pvec(2) * pzpt, nsvahl * pvec(1) / pt * sqh );
      }
      else
      {
        vc[3] = cxmake( -hel * sqh, 0. );
            // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
            //vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
        vc[4] = cxmake( 0., nsvahl * (pvec(3) < 0 ? -sqh : sqh ) );
      }
    }
  }
  else
  {
    const fptype_sv& pp = pvec(0); // NB: rewrite the following as in Fortran, using pp instead of pvec0
    const fptype_sv pt = fpsqrt( ( pvec(1) * pvec(1) ) + ( pvec(2) * pvec(2) ) );
    vc[2] = cxzero_sv();
    vc[5] = cxmake( hel * pt / pp * sqh, 0. );
    if ( pt != 0. )
    {
      const fptype pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
      vc[3] = cxmake( -pvec(1) * pzpt, -nsv * pvec(2) / pt * sqh );
      vc[4] = cxmake( -pvec(2) * pzpt, nsv * pvec(1) / pt * sqh );
    }
    else
    {
      vc[3] = cxmake( -hel * sqh, 0. );
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          //vc[4] = cxmake( 0, nsv * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV why abs here?
      vc[4] = cxmake( 0., nsv * ( pvec(3) < 0. ? -sqh : sqh ) ); // AV: removed an abs here
    }
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION  void sxxxxx(const T& pvec,
                             const fptype& smass,
                             const int& nhel,
                             const int& nss,
                             cxtype_sv sc[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  sc[2] = cxmake( 1 + fptype_sv{0}, 0 );
      sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
      sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void oxxxxx(const T& pvec,
                            const fptype fmass,
                            const int nhel,
                            const int nsf,
                            cxtype_sv fo[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
  fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
  const int nh = nhel * nsf;
  if ( fmass != 0. )
  {
    const fptype_sv pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
    if ( pp == 0. )
    {
      // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
      fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
      //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
      sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
      const int ip = -( ( 1 - nh ) / 2 ) * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
      const int im = ( 1 + nh ) / 2 * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
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
  }
  else
  {
    const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. ) and ( pvec2 == 0. ) and ( pvec3 < 0. ),
                                        0, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
    const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ),
                                cxternary( ( sqp0p3 == 0. ),
                                          cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 ),
                                          cxmake( (fptype)nh * pvec1, -pvec2 ) / sqp0p3 ) };
    if ( nh == 1 )
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void opzxxx(const T& pvec,
                             const int& nhel,                  // input: -1 or +1 (helicity of fermion)
                             const int& nsf,                  // input: +1 (particle) or -1 (antiparticle)
                             cxtype_sv fo[])
{

  const fptype pvec3 = pvec(3);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void omzxxx(const T& pvec,
                             const int& nhel,
                             const int& nsf,
                             cxtype_sv fo[])
{

  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( -pvec3 * (fptype)nsf, pvec3 * (fptype)nsf ); // remember pvec0 == -pvec3
  fo[1] = cxzero_sv();
  const int nh = nhel * nsf;
  const cxtype_sv chi1 = cxmake( -nhel, 0. ) * fpsqrt( -2. * pvec3 );

  if ( nh == 1 )
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void oxzxxx(const T& pvec,
                             const int& nhel,
                             const int& nsf,
                             cxtype_sv fo[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
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
  return;
}



  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVV1_0( const cxtype V1[],
               const cxtype V2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVV1P0_1( const cxtype V2[],
                 const cxtype V3[],
                 const cxtype COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype V1[] );

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_1( const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               const fptype M1,
               const fptype W1,
               cxtype F1[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_2( const cxtype F1[],
               const cxtype V3[],
               const cxtype COUP,
               const fptype M2,
               const fptype W2,
               cxtype F2[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1P0_3( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV1P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV3P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV4P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] );

  //--------------------------------------------------------------------------


// --------------------------------------------------------------------------



  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVV1_0( const cxtype V1[],
               const cxtype V2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const fptype_sv P1[4] = { +cxreal( V1[0] ), +cxreal( V1[1] ), +cximag( V1[1] ), +cximag( V1[0] ) };
    const fptype_sv P2[4] = { +cxreal( V2[0] ), +cxreal( V2[1] ), +cximag( V2[1] ), +cximag( V2[0] ) };
    const fptype_sv P3[4] = { +cxreal( V3[0] ), +cxreal( V3[1] ), +cximag( V3[1] ), +cximag( V3[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP1 = ( V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5] );
    const cxtype_sv TMP2 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP3 = ( V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5] );
    const cxtype_sv TMP4 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP5 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv TMP7 = ( V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3] );
    const cxtype_sv TMP8 = ( V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3] );
    (*vertex) = COUP * ( TMP1 * ( - cI * ( TMP0 ) + cI * ( TMP2 ) ) + ( TMP3 * ( +cI * ( TMP4 )- cI * ( TMP5 ) ) + TMP6 * ( - cI * ( TMP7 ) + cI * ( TMP8 ) ) ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVV1P0_1( const cxtype V2[],
                 const cxtype V3[],
                 const cxtype COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype V1[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const fptype_sv P2[4] = { +cxreal( V2[0] ), +cxreal( V2[1] ), +cximag( V2[1] ), +cximag( V2[0] ) };
    const fptype_sv P3[4] = { +cxreal( V3[0] ), +cxreal( V3[1] ), +cximag( V3[1] ), +cximag( V3[0] ) };
    V1[0] = + V2[0] + V3[0];
    V1[1] = + V2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP2 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP4 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP5 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( TMP6 * ( - cI * ( P2[0] ) + cI * ( P3[0] ) ) + ( V2[2] * ( - cI * ( TMP0 ) + cI * ( TMP2 ) ) + V3[2] * ( +cI * ( TMP4 )- cI * ( TMP5 ) ) ) );
    V1[3] = denom * ( TMP6 * ( - cI * ( P2[1] ) + cI * ( P3[1] ) ) + ( V2[3] * ( - cI * ( TMP0 ) + cI * ( TMP2 ) ) + V3[3] * ( +cI * ( TMP4 )- cI * ( TMP5 ) ) ) );
    V1[4] = denom * ( TMP6 * ( - cI * ( P2[2] ) + cI * ( P3[2] ) ) + ( V2[4] * ( - cI * ( TMP0 ) + cI * ( TMP2 ) ) + V3[4] * ( +cI * ( TMP4 )- cI * ( TMP5 ) ) ) );
    V1[5] = denom * ( TMP6 * ( - cI * ( P2[3] ) + cI * ( P3[3] ) ) + ( V2[5] * ( - cI * ( TMP0 ) + cI * ( TMP2 ) ) + V3[5] * ( +cI * ( TMP4 )- cI * ( TMP5 ) ) ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const cxtype_sv TMP9 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + ( F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    (*vertex) = COUP * - cI * TMP9;
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_1( const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               const fptype M1,
               const fptype W1,
               cxtype F1[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    F1[0] = + F2[0] + V3[0];
    F1[1] = + F2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( F1[0] ), -cxreal( F1[1] ), -cximag( F1[1] ), -cximag( F1[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    F1[2] = denom * cI * ( F2[2] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3]- cI * ( V3[4] ) ) + ( P1[2] * ( +cI * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[1] * (- one) * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + M1 * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) ) );
    F1[3] = denom * (- cI) * ( F2[2] * ( P1[0] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( V3[2] - V3[5] ) + ( P1[2] * ( - cI * ( V3[2] ) + cI * ( V3[5] ) ) + P1[3] * ( V3[3]- cI * ( V3[4] ) ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P1[2] * ( +cI * ( V3[3] ) - V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + M1 * ( F2[4] * ( -V3[3] + cI * ( V3[4] ) ) + F2[5] * ( -V3[2] + V3[5] ) ) ) );
    F1[4] = denom * (- cI) * ( F2[4] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[2] * (- one) * ( +cI * ( V3[3] ) + V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + ( F2[5] * ( P1[0] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( -V3[2] + V3[5] ) + ( P1[2] * ( - cI * ( V3[2] ) + cI * ( V3[5] ) ) - P1[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + M1 * ( F2[2] * ( -V3[2] + V3[5] ) + F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) );
    F1[5] = denom * cI * ( F2[4] * ( P1[0] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F2[5] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[2] * ( - cI * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + M1 * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_2( const cxtype F1[],
               const cxtype V3[],
               const cxtype COUP,
               const fptype M2,
               const fptype W2,
               cxtype F2[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    F2[0] = + F1[0] + V3[0];
    F2[1] = + F1[1] + V3[1];
    const fptype_sv P2[4] = { -cxreal( F2[0] ), -cxreal( F2[1] ), -cximag( F2[1] ), -cximag( F2[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P2[0] * P2[0] ) - ( P2[1] * P2[1] ) - ( P2[2] * P2[2] ) - ( P2[3] * P2[3] ) - M2 * ( M2 - cI * W2 ) );
    F2[2] = denom * cI * ( F1[2] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[2] * ( +cI * ( V3[3] ) - V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + ( F1[3] * ( P2[0] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[1] * ( -V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] )- cI * ( V3[5] ) ) + P2[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) ) + M2 * ( F1[4] * ( V3[2] - V3[5] ) + F1[5] * ( -V3[3] + cI * ( V3[4] ) ) ) ) );
    F2[3] = denom * (- cI) * ( F1[2] * ( P2[0] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[1] * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) - P2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F1[3] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[2] * ( +cI * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + M2 * ( F1[4] * ( V3[3] + cI * ( V3[4] ) ) - F1[5] * ( V3[2] + V3[5] ) ) ) );
    F2[4] = denom * (- cI) * ( F1[4] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3] + cI * ( V3[4] ) ) + ( P2[2] * ( - cI * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[1] * (- one) * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) + P2[3] * ( V3[3]- cI * ( V3[4] ) ) ) ) ) + M2 * ( F1[2] * (- one) * ( V3[2] + V3[5] ) + F1[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) );
    F2[5] = denom * cI * ( F1[4] * ( P2[0] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[1] * ( V3[2] - V3[5] ) + ( P2[2] * ( +cI * ( V3[2] )- cI * ( V3[5] ) ) + P2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * ( -V3[3] + cI * ( V3[4] ) ) + ( P2[2] * (- one) * ( +cI * ( V3[3] ) + V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + M2 * ( F1[2] * ( V3[3] + cI * ( V3[4] ) ) + F1[3] * ( V3[2] - V3[5] ) ) ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1P0_3( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype_sv denom = COUP / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * (- cI) * ( F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3] );
    V3[3] = denom * (- cI) * ( -F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] );
    V3[4] = denom * (- cI) * ( - cI * ( F1[2] * F2[5] + F1[5] * F2[2] ) + cI * ( F1[3] * F2[4] + F1[4] * F2[3] ) );
    V3[5] = denom * (- cI) * ( -F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2] );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV1P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP10 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( - cI * ( TMP6 * V4[2] ) + cI * ( V3[2] * TMP10 ) );
    V1[3] = denom * ( - cI * ( TMP6 * V4[3] ) + cI * ( V3[3] * TMP10 ) );
    V1[4] = denom * ( - cI * ( TMP6 * V4[4] ) + cI * ( V3[4] * TMP10 ) );
    V1[5] = denom * ( - cI * ( TMP6 * V4[5] ) + cI * ( V3[5] * TMP10 ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV3P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP11 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( - cI * ( TMP6 * V4[2] ) + cI * ( V2[2] * TMP11 ) );
    V1[3] = denom * ( - cI * ( TMP6 * V4[3] ) + cI * ( V2[3] * TMP11 ) );
    V1[4] = denom * ( - cI * ( TMP6 * V4[4] ) + cI * ( V2[4] * TMP11 ) );
    V1[5] = denom * ( - cI * ( TMP6 * V4[5] ) + cI * ( V2[5] * TMP11 ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  
  KOKKOS_INLINE_FUNCTION
  void VVVV4P0_1( const cxtype V2[],
                  const cxtype V3[],
                  const cxtype V4[],
                  const cxtype COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype V1[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP10 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP11 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( - cI * ( V3[2] * TMP10 ) + cI * ( V2[2] * TMP11 ) );
    V1[3] = denom * ( - cI * ( V3[3] * TMP10 ) + cI * ( V2[3] * TMP11 ) );
    V1[4] = denom * ( - cI * ( V3[4] * TMP10 ) + cI * ( V2[4] * TMP11 ) );
    V1[5] = denom * ( - cI * ( V3[5] * TMP10 ) + cI * ( V2[5] * TMP11 ) );
    return;
  }

  //--------------------------------------------------------------------------



} // end namespace MG5_sm

#endif // HelAmps_sm_H
