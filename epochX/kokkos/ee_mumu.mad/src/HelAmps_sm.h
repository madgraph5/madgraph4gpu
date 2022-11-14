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

namespace MG5_sm {



//--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------
  // Memory function not needed in Kokkos

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ixxxxx(const mom_t& pvec,
          const fptype fmass,             // input: fermion mass
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fi);

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ipzxxx(const mom_t& pvec,
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fi);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void imzxxx(const mom_t& pvec, 
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fi);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ixzxxx(const mom_t& pvec,
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fi);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void vxxxxx(const mom_t& pvec,
          const fptype vmass,             // input: vector boson mass
          const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
          const int nsv,                  // input: +1 (final) or -1 (initial)
          cxtype* vc);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void sxxxxx(const mom_t& pvec,
          //const fptype,                 // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
          //const int,                    // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
          const int nss,                  // input: +1 (final) or -1 (initial)
          cxtype* sc);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void oxxxxx(const mom_t& pvec,
          const fptype fmass,             // input: fermion mass
          const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fo);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void opzxxx(const mom_t& pvec,
          //const fptype fmass,           // ASSUME fermion mass==0
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fo);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void omzxxx(const mom_t& pvec,
          //const fptype fmass,           // ASSUME fermion mass==0
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fo);

  //--------------------------------------------------------------------------
  
  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void oxzxxx(const mom_t& pvec,
          //const fptype fmass,           // ASSUME fermion mass==0
          const int nhel,                 // input: -1 or +1 (helicity of fermion)
          const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
          cxtype* fo);

  //==========================================================================

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ixxxxx(const mom_t& pvec,
                            const fptype fmass,             // input: fermion mass
                            const int nhel,                 // input: -1 or +1 (helicity of fermion)
                            const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
                            cxtype* fi)
{
  const fptype& pvec0 = pvec(0);
  const fptype& pvec1 = pvec(1);
  const fptype& pvec2 = pvec(2);
  const fptype& pvec3 = pvec(3);
  fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf );
  fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf );
  const int nh = nhel * nsf;
  if ( fmass != 0. )
  {
      const fptype pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
    if ( pp == 0. )
    {
      // NB: Do not use "fpabs" for floats! It returns an integer with no build warning! Use fpabs!
      fptype sqm[2] = { fpsqrt( fpabs( fmass ) ), 0. }; // possibility of negative fermion masses
          //sqm[1] = ( fmass < 0. ? -fpabs( sqm[0] ) : fpabs( sqm[0] ) ); // AV: why fpabs here?
      sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an fpabs here
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
  }
  else
  { 
      const fptype sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                          fptype{0}, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
    const cxtype chi[2] = { cxmake( sqp0p3, 0. ), cxternary( ( sqp0p3 == 0. ),
                                                                  cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0 ), 0. ),
                                                                  cxmake( (fptype)nh * pvec1, pvec2 ) / sqp0p3 ) };
    if ( nh == 1 )
    {
      fi[2] = cxtype( 0, 0 );
      fi[3] = cxtype( 0, 0 );
      fi[4] = chi[0];
      fi[5] = chi[1];
    }
    else
    {
      fi[2] = chi[1];
      fi[3] = chi[0];
      fi[4] = cxtype( 0, 0 );
      fi[5] = cxtype( 0, 0 );
    }
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ipzxxx(const mom_t& pvec,
                            const int nhel,
                            const int nsf,
                            cxtype* fi)
{
  const fptype& pvec3 = pvec(3);
  
  fi[0] = cxmake(-pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf);
  fi[1] = cxtype( 0, 0 );
  const int nh = nhel * nsf;
  const cxtype sqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
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
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void imzxxx(const mom_t& pvec, 
                            const int nhel,
                            const int nsf,
                            cxtype* fi)
{

  const fptype pvec3 = pvec(3);
  fi[0] = cxmake( pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
  fi[1] = cxtype( 0, 0 );
  const int nh = nhel * nsf;
  const cxtype chi = cxmake( -(fptype)nhel * fpsqrt( -2. * pvec3 ), 0. );
  fi[3] = cxtype( 0, 0 );
  fi[4] = cxtype( 0, 0 );
  if ( nh == 1 )
  {
    fi[2] = cxtype( 0, 0 );
     fi[5] = chi;
  }
  else
  {
     fi[2] = chi;
     fi[5] = cxtype( 0, 0 );
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void ixzxxx(const mom_t& pvec,
                            const int nhel,
                            const int nsf,
                            cxtype* fi)
{
  const fptype& pvec0 = pvec(0);
  const fptype& pvec1 = pvec(1);
  const fptype& pvec2 = pvec(2);
  const fptype& pvec3 = pvec(3);

  //fi[0] = cxmake( -pvec0 * nsf, -pvec2 * nsf ); // AV: BUG! not the same as ixxxxx
  //fi[1] = cxmake( -pvec0 * nsf, -pvec1 * nsf ); // AV: BUG! not the same as ixxxxx
  fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf ); // AV: BUG FIX
  fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf ); // AV: BUG FIX
  const int nh = nhel * nsf;
  //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
  const fptype sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
  const cxtype chi0 = cxmake( sqp0p3, 0. );
  const cxtype chi1 = cxmake( (fptype)nh * pvec1/sqp0p3, pvec2/sqp0p3 );
  if ( nh == 1 )
  {
    fi[2] = cxtype( 0, 0 );
    fi[3] = cxtype( 0, 0 );
    fi[4] = chi0;
    fi[5] = chi1;
  }
  else
  {
    fi[2] = chi1;
    fi[3] = chi0;
    fi[4] = cxtype( 0, 0 );
    fi[5] = cxtype( 0, 0 );
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void vxxxxx(const mom_t& pvec,
                            const fptype vmass,
                            const int nhel,
                            const int nsv,
                            cxtype* vc) 
{
  
  const fptype& pvec0 = pvec(0);
  const fptype& pvec1 = pvec(1);
  const fptype& pvec2 = pvec(2);
  const fptype& pvec3 = pvec(3);
  const fptype sqh = fpsqrt( 0.5 );
  const fptype hel = nhel;
  const int nsvahl = nsv * fpabs(hel);
  vc[0] = cxmake( pvec0 * (fptype)nsv, pvec3 * (fptype)nsv );
  vc[1] = cxmake( pvec1 * (fptype)nsv, pvec2 * (fptype)nsv );
  if ( vmass != 0. )
  {
    const int nsvahl = nsv * fpabs( hel );
    const fptype pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
    const fptype pp = fpmin( pvec0, fpsqrt( pt2 + ( pvec3 * pvec3 ) ) );
    const fptype pt = fpmin( pp, fpsqrt( pt2 ) );
    const fptype hel0 = 1. - fpabs( hel );
    if ( pp == 0. )
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
      if ( pt != 0. )
      {
          const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
          vc[3] = cxmake( hel0 * pvec1 * emp - pvec1 * pzpt, -nsvahl * pvec2 / pt * sqh );
          vc[4] = cxmake( hel0 * pvec2 * emp - pvec2 * pzpt, nsvahl * pvec1 / pt * sqh );
      }
      else
      {
        vc[3] = cxmake( -hel * sqh, 0. );
          vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
      }
    }
  }
  else
  {
      const fptype& pp = pvec0; // NB: rewrite the following as in Fortran, using pp instead of pvec0
      const fptype pt = fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) );
    vc[2] = cxtype( 0, 0 );
    vc[5] = cxmake( hel * pt / pp * sqh, 0. );
    if ( pt != 0. )
    {
        const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
        vc[3] = cxmake( -pvec1 * pzpt, -nsv * pvec2 / pt * sqh );
        vc[4] = cxmake( -pvec2 * pzpt, nsv * pvec1 / pt * sqh );
    }
    else
    {
      vc[3] = cxmake( -hel * sqh, 0. );
        vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
    }
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void sxxxxx(const mom_t& pvec,
                             const fptype& smass,
                             const int nhel,
                             const int nss,
                             cxtype* sc)
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  sc[2] = cxmake( 1 + fptype{0}, 0 );
      sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
      sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
template<typename mom_t>
KOKKOS_INLINE_FUNCTION void oxxxxx(const mom_t& pvec,
                            const fptype fmass,
                            const int nhel,
                            const int nsf,
                            cxtype* fo)
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
    const fptype pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
    if ( pp == 0. )
    {
      fptype sqm[2] = { fpsqrt( fpabs( fmass ) ), 0. }; // possibility of negative fermion masses
      //sqm[1] = ( fmass < 0. ? -fpabs( sqm[0] ) : fpabs( sqm[0] ) ); // AV: why fpabs here?
      sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an fpabs here
      const int ip = -( ( 1 - nh ) / 2 ) * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
      const int im = ( 1 + nh ) / 2 * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
      fo[2] = cxmake( im * sqm[iabs( ip )], 0 );
      fo[3] = cxmake( ip * nsf * sqm[iabs( ip )], 0 );
      fo[4] = cxmake( im * nsf * sqm[iabs( im )], 0 );
      fo[5] = cxmake( ip * sqm[iabs( im )], 0 );
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
    const fptype sqp0p3 = fpternary( ( pvec1 == 0. ) and ( pvec2 == 0. ) and ( pvec3 < 0. ),
                                        0, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
    const cxtype chi[2] = { cxmake( sqp0p3, 0. ),
                                cxternary( ( sqp0p3 == 0. ),
                                          cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 ),
                                          cxmake( (fptype)nh * pvec1, -pvec2 ) / sqp0p3 ) };
    if ( nh == 1 )
    {
      fo[2] = chi[0];
      fo[3] = chi[1];
      fo[4] = cxtype( 0, 0 );
      fo[5] = cxtype( 0, 0 );
    }
    else
    {
      fo[2] = cxtype( 0, 0 );
      fo[3] = cxtype( 0, 0 );
      fo[4] = chi[1];
      fo[5] = chi[0];
    }
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void opzxxx(const mom_t& pvec,
                             const int nhel,                  // input: -1 or +1 (helicity of fermion)
                             const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
                             cxtype* fo)
{

  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( pvec3 * (fptype)nsf, pvec3 * (fptype)nsf );
  fo[1] = cxtype( 0, 0 );
  const int nh = nhel * nsf;
  const cxtype csqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
  fo[3] = cxtype( 0, 0 );
  fo[4] = cxtype( 0, 0 );
  if ( nh == 1 )
  {
    fo[2] = csqp0p3;
    fo[5] = cxtype( 0, 0 );
  }
  else
  {
    fo[2] = cxtype( 0, 0 );
    fo[5] = csqp0p3;
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void omzxxx(const mom_t& pvec,
                             const int nhel,
                             const int nsf,
                             cxtype* fo)
{

  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( -pvec3 * (fptype)nsf, pvec3 * (fptype)nsf ); // remember pvec0 == -pvec3
  fo[1] = cxtype( 0, 0 );
  const int nh = nhel * nsf;
  const cxtype chi1 = cxmake( -nhel, 0. ) * fpsqrt( -2. * pvec3 );
  if ( nh == 1 )
  {
    fo[2] = cxtype( 0, 0 );
    fo[3] = chi1;
    fo[4] = cxtype( 0, 0 );
    fo[5] = cxtype( 0, 0 );
  }
  else
  {
    fo[2] = cxtype( 0, 0 );
    fo[3] = cxtype( 0, 0 );
    fo[4] = chi1;
    //fo[5] = chi1; // AV: BUG!
    fo[5] = cxtype( 0, 0 ); // AV: BUG FIX
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename mom_t>
KOKKOS_INLINE_FUNCTION  void oxzxxx(const mom_t& pvec,
                             const int nhel,
                             const int nsf,
                             cxtype* fo)
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
  fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
  const int nh = nhel * nsf;
  //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
  const fptype sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
  const cxtype chi0 = cxmake( sqp0p3, 0. );
  const cxtype chi1 = cxmake( (fptype)nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
  if ( nh == 1 )
  {
    fo[2] = chi0;
    fo[3] = chi1;
    fo[4] = cxtype( 0, 0 );
    fo[5] = cxtype( 0, 0 );
  }
  else
  {
    fo[2] = cxtype( 0, 0 );
    fo[3] = cxtype( 0, 0 );
    fo[4] = chi1;
    fo[5] = chi0;
  }
  return;
}



  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] );

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

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] );

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV4_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV4_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] );

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_4_0( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype V3[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype vertex[] );

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_4_3( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] );

  //--------------------------------------------------------------------------


// --------------------------------------------------------------------------



  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV1_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const cxtype TMP0 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) + ( F1[3] * ( F2[4] * ( V3[3] - cI * V3[4] ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * V3[4] ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * V3[4] ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    ( *vertex ) = COUP * -cI * TMP0;
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
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype denom = COUP / ( ( P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * ( -cI ) * ( F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3] );
    V3[3] = denom * ( -cI ) * ( -F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] );
    V3[4] = denom * ( -cI ) * ( -cI * ( F1[2] * F2[5] + F1[5] * F2[2] ) + cI * ( F1[3] * F2[4] + F1[4] * F2[3] ) );
    V3[5] = denom * ( -cI ) * ( -F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2] );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const cxtype TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) + F1[3] * ( F2[4] * ( V3[3] - cI * V3[4] ) + F2[5] * ( V3[2] - V3[5] ) ) );
    ( *vertex ) = COUP * -cI * TMP1;
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * P3[2] ) ) + F1[3] * ( F2[4] * ( P3[1] - cI * P3[2] ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype denom = COUP / ( ( P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * ( -cI ) * ( F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2 );
    V3[3] = denom * ( -cI ) * ( -F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP2 );
    V3[4] = denom * ( -cI ) * ( -cI * ( F1[2] * F2[5] ) + cI * ( F1[3] * F2[4] ) - P3[2] * OM3 * TMP2 );
    V3[5] = denom * ( -cI ) * ( -F1[2] * F2[4] - P3[3] * OM3 * TMP2 + F1[3] * F2[5] );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV4_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    const cxtype TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) + F1[3] * ( F2[4] * ( V3[3] - cI * V3[4] ) + F2[5] * ( V3[2] - V3[5] ) ) );
    const cxtype TMP3 = ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * V3[4] ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * V3[4] ) + F2[3] * ( V3[2] + V3[5] ) ) );
    ( *vertex ) = COUP * ( -one ) * ( +cI * TMP1 + ( two * cI ) * TMP3 );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV4_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const cxtype TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * P3[2] ) ) + F1[3] * ( F2[4] * ( P3[1] - cI * P3[2] ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype TMP4 = ( F1[4] * ( F2[2] * ( P3[0] - P3[3] ) - F2[3] * ( P3[1] + cI * P3[2] ) ) + F1[5] * ( F2[2] * ( -P3[1] + cI * P3[2] ) + F2[3] * ( P3[0] + P3[3] ) ) );
    const cxtype denom = COUP / ( ( P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * ( -two * cI ) * ( OM3 * -half * P3[0] * ( TMP2 + two * TMP4 ) + ( +half * ( F1[2] * F2[4] + F1[3] * F2[5] ) + F1[4] * F2[2] + F1[5] * F2[3] ) );
    V3[3] = denom * ( -two * cI ) * ( OM3 * -half * P3[1] * ( TMP2 + two * TMP4 ) + ( -half * ( F1[2] * F2[5] + F1[3] * F2[4] ) + F1[4] * F2[3] + F1[5] * F2[2] ) );
    V3[4] = denom * ( two * cI ) * ( OM3 * half * P3[2] * ( TMP2 + two * TMP4 ) + ( half * cI * ( F1[2] * F2[5] ) - half * cI * ( F1[3] * F2[4] ) - cI * ( F1[4] * F2[3] ) + cI * ( F1[5] * F2[2] ) ) );
    V3[5] = denom * ( two * cI ) * ( OM3 * half * P3[3] * ( TMP2 + two * TMP4 ) + ( +half * ( F1[2] * F2[4] ) - half * ( F1[3] * F2[5] ) - F1[4] * F2[2] + F1[5] * F2[3] ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_4_0( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype V3[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype vertex[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    const cxtype TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * V3[4] ) ) + F1[3] * ( F2[4] * ( V3[3] - cI * V3[4] ) + F2[5] * ( V3[2] - V3[5] ) ) );
    const cxtype TMP3 = ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * V3[4] ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * V3[4] ) + F2[3] * ( V3[2] + V3[5] ) ) );
    ( *vertex ) = ( -one ) * ( COUP2 * ( +cI * TMP1 + ( two * cI ) * TMP3 ) + cI * ( TMP1 * COUP1 ) );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  
  KOKKOS_INLINE_FUNCTION
  void FFV2_4_3( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] )
  {
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const cxtype TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * P3[2] ) ) + F1[3] * ( F2[4] * ( P3[1] - cI * P3[2] ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype TMP4 = ( F1[4] * ( F2[2] * ( P3[0] - P3[3] ) - F2[3] * ( P3[1] + cI * P3[2] ) ) + F1[5] * ( F2[2] * ( -P3[1] + cI * P3[2] ) + F2[3] * ( P3[0] + P3[3] ) ) );
    const cxtype denom = one / ( ( P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * ( -two * cI ) * ( COUP2 * ( OM3 * -half * P3[0] * ( TMP2 + two * TMP4 ) + ( +half * ( F1[2] * F2[4] + F1[3] * F2[5] ) + F1[4] * F2[2] + F1[5] * F2[3] ) ) + half * ( COUP1 * ( F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2 ) ) );
    V3[3] = denom * ( -two * cI ) * ( COUP2 * ( OM3 * -half * P3[1] * ( TMP2 + two * TMP4 ) + ( -half * ( F1[2] * F2[5] + F1[3] * F2[4] ) + F1[4] * F2[3] + F1[5] * F2[2] ) ) - half * ( COUP1 * ( F1[2] * F2[5] + F1[3] * F2[4] + P3[1] * OM3 * TMP2 ) ) );
    V3[4] = denom * cI * ( COUP2 * ( OM3 * P3[2] * ( TMP2 + two * TMP4 ) + ( +cI * ( F1[2] * F2[5] ) - cI * ( F1[3] * F2[4] ) + ( -two * cI ) * ( F1[4] * F2[3] ) + ( two * cI ) * ( F1[5] * F2[2] ) ) ) + COUP1 * ( +cI * ( F1[2] * F2[5] ) - cI * ( F1[3] * F2[4] ) + P3[2] * OM3 * TMP2 ) );
    V3[5] = denom * ( two * cI ) * ( COUP2 * ( OM3 * half * P3[3] * ( TMP2 + two * TMP4 ) + ( +half * ( F1[2] * F2[4] ) - half * ( F1[3] * F2[5] ) - F1[4] * F2[2] + F1[5] * F2[3] ) ) + half * ( COUP1 * ( F1[2] * F2[4] + P3[3] * OM3 * TMP2 - F1[3] * F2[5] ) ) );
    return;
  }

  //--------------------------------------------------------------------------



} // end namespace MG5_sm

#endif // HelAmps_sm_H
