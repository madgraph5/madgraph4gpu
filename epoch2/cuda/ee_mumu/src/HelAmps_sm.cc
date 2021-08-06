//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

mgDebugDeclare();

namespace MG5_sm
{
  __device__
  inline const fptype& pIparIp4Ievt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    //return allmomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
    return momenta[ipagM][ipar][ip4][ieppM];
  }

  //--------------------------------------------------------------------------

  __global__
  void ixxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly slower, in c++
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      //const fptype pvec0 = fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ); // AV: BUG?! (NOT AS IN THE FORTRAN)
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt ); // AV: BUG FIX (DO AS IN THE FORTRAN)
      fi[0] = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi[1] = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      if ( fmass != 0. )
      {
        const fptype pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
        if ( pp == 0. )
        {
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
          //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
          sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
          const int ip = ( 1 + nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          const int im = ( 1 - nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          fi[2] = ip * sqm[ip];
          fi[3] = im * nsf * sqm[ip];
          fi[4] = ip * nsf * sqm[im];
          fi[5] = im * sqm[im];
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
        const fptype sqp0p3 = ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ?
                                0. : fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * nsf );
        const cxtype chi[2] = { cxmake( sqp0p3, 0. ), ( sqp0p3 == 0. ?
                                                        cxmake( -nhel * fpsqrt( 2. * pvec0 ), 0. ) :
                                                        cxmake( nh * pvec1, pvec2 ) / sqp0p3 ) };
        if ( nh == 1 )
        {
          fi[2] = cxmake( 0., 0. );
          fi[3] = cxmake( 0., 0. );
          fi[4] = chi[0];
          fi[5] = chi[1];
        }
        else
        {
          fi[2] = chi[1];
          fi[3] = chi[0];
          fi[4] = cxmake( 0., 0. );
          fi[5] = cxmake( 0., 0. );
        }
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void ipzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fi[0] = cxmake( -pvec3 * nsf, -pvec3 * nsf );
      fi[1] = cxmake( 0., 0. );
      const int nh = nhel * nsf;
      const cxtype sqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * nsf, 0. );
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void imzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fi[0] = cxmake( pvec3 * nsf, -pvec3 * nsf );
      fi[1] = cxmake( 0., 0. );
      const int nh = nhel * nsf;
      const cxtype chi = cxmake( -nhel * fpsqrt( -2. * pvec3 ), 0. );
      fi[3] = cxmake( 0., 0. );
      fi[4] = cxmake( 0., 0. );
      if ( nh == 1 )
      {
        fi[2] = cxmake( 0., 0. );
        fi[5] = chi;
      }
      else
      {
        fi[2] = chi;
        fi[5] = cxmake( 0., 0. );
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void ixzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      //fi[0] = cxmake( -pvec0 * nsf, -pvec2 * nsf ); // AV: BUG! not the same as ixxxxx
      //fi[1] = cxmake( -pvec0 * nsf, -pvec1 * nsf ); // AV: BUG! not the same as ixxxxx
      fi[0] = cxmake( -pvec0 * nsf, -pvec3 * nsf ); // AV: BUG FIX
      fi[1] = cxmake( -pvec1 * nsf, -pvec2 * nsf ); // AV: BUG FIX
      const int nh = nhel * nsf;
      //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
      const fptype sqp0p3 = fpsqrt( pvec0 + pvec3 ) * nsf;
      const cxtype chi0 = cxmake( sqp0p3, 0. );
      const cxtype chi1 = cxmake( nh * pvec1/sqp0p3, pvec2/sqp0p3 );
      if ( nh == 1 )
      {
        fi[2] = cxmake( 0., 0. );
        fi[3] = cxmake( 0., 0. );
        fi[4] = chi0;
        fi[5] = chi1;
      }
      else
      {
        fi[2] = chi1;
        fi[3] = chi0;
        fi[4] = cxmake( 0., 0. );
        fi[5] = cxmake( 0., 0. );
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void vxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,           // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,            // input: +1 (final) or -1 (initial)
               cxtype* vc,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      const fptype sqh = fpsqrt( 0.5 );
      const fptype hel = nhel;
      vc[0] = cxmake( pvec0 * nsv, pvec3 * nsv );
      vc[1] = cxmake( pvec1 * nsv, pvec2 * nsv );
      if ( vmass != 0. )
      {
        const int nsvahl = nsv * std::abs( hel );
        const fptype pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
        const fptype pp = fpmin( pvec0, fpsqrt( pt2 + ( pvec3 * pvec3 ) ) );
        const fptype pt = fpmin( pp, fpsqrt( pt2 ) );
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
            // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
            //vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
            vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
          }
        }
      }
      else
      {
        const fptype& pp = pvec0; // NB: rewrite the following  as in Fortran, using pp instead of pvec0
        const fptype pt = fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) );
        vc[2] = cxmake( 0., 0. );
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
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          //vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
          vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
        }
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void sxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype,             // WARNING: "smass" unused (missing in Fortran)
               const int,                // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,            // input: +1 (final) or -1 (initial)
               cxtype sc[3],             // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      sc[2] = cxmake( 1., 0. );
      sc[0] = cxmake( pvec0 * nss, pvec3 * nss );
      sc[1] = cxmake( pvec1 * nss, pvec2 * nss );
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void oxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly faster, in c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fo[0] = cxmake( pvec0 * nsf, pvec3 * nsf );
      fo[1] = cxmake( pvec1 * nsf, pvec2 * nsf );
      const int nh = nhel * nsf;
      if ( fmass != 0. )
      {
        const fptype pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
        if ( pp == 0. )
        {
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
          //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
          sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
          const int ip = -( ( 1 - nh ) / 2 ) * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          const int im = ( 1 + nh ) / 2 * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          fo[2] = im * sqm[std::abs( ip )];
          fo[3] = ip * nsf * sqm[std::abs( ip )];
          fo[4] = im * nsf * sqm[std::abs( im )];
          fo[5] = ip * sqm[std::abs( im )];
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
        const fptype sqp0p3 = ( ( pvec1 == 0. ) and ( pvec2 == 0. ) and ( pvec3 < 0. )
                                ? 0. : fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * nsf );
        const cxtype chi[2] = { cxmake( sqp0p3, 0. ),
                                ( ( sqp0p3 == 0. ) ? cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 )
                                  : cxmake( nh * pvec1, -pvec2 ) / sqp0p3 ) };
        if ( nh == 1 )
        {
          fo[2] = chi[0];
          fo[3] = chi[1];
          fo[4] = cxmake( 0., 0. );
          fo[5] = cxmake( 0., 0. );
        }
        else
        {
          fo[2] = cxmake( 0., 0. );
          fo[3] = cxmake( 0., 0. );
          fo[4] = chi[1];
          fo[5] = chi[0];
        }
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void opzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fo[0] = cxmake( pvec3 * nsf, pvec3 * nsf );
      fo[1] = cxmake( 0., 0. );
      const int nh = nhel * nsf;
      const cxtype csqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * nsf, 0. );
      fo[3] = cxmake( 0., 0. );
      fo[4] = cxmake( 0., 0. );
      if ( nh == 1 )
      {
        fo[2] = csqp0p3;
        fo[5] = cxmake( 0., 0. );
      }
      else
      {
        fo[2] = cxmake( 0., 0. );
        fo[5] = csqp0p3;
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void omzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
#endif
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fo[0] = cxmake( -pvec3 * nsf, pvec3 * nsf ); // remember pvec0 == -pvec3
      fo[1] = cxmake( 0., 0. );
      const int nh = nhel * nsf;
      const cxtype chi1 = cxmake( -nhel, 0. ) * fpsqrt( -2. * pvec3 );
      if ( nh == 1 )
      {
        fo[2] = cxmake( 0., 0. );
        fo[3] = chi1;
        fo[4] = cxmake( 0., 0. );
        fo[5] = cxmake( 0., 0. );
      }
      else
      {
        fo[2] = cxmake( 0., 0. );
        fo[3] = cxmake( 0., 0. );
        fo[4] = chi1;
        //fo[5] = chi1; // AV: BUG!
        fo[5] = cxmake( 0., 0. ); // AV: BUG FIX
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void oxzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,           // input: -1 or +1 (helicity of fermion)
               const int nsf,            // input: +1 (particle) or -1 (antiparticle)
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      fo[0] = cxmake( pvec0 * nsf, pvec3 * nsf );
      fo[1] = cxmake( pvec1 * nsf, pvec2 * nsf );
      const int nh = nhel * nsf;
      //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
      const fptype sqp0p3 = fpsqrt( pvec0 + pvec3 ) * nsf;
      const cxtype chi0 = cxmake( sqp0p3, 0. );
      const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
      if ( nh == 1 )
      {
        fo[2] = chi0;
        fo[3] = chi1;
        fo[4] = cxmake( 0., 0. );
        fo[5] = cxmake( 0., 0. );
      }
      else
      {
        fo[2] = cxmake( 0., 0. );
        fo[3] = cxmake( 0., 0. );
        fo[4] = chi1;
        fo[5] = chi0;
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV1_0( const cxtype F1[],    // input: wavefunction1[6]
               const cxtype F2[],    // input: wavefunction2[6]
               const cxtype V3[],    // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )      // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    const cxtype TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
                         (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
                          (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
                           F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])))));
    (*vertex) = COUP * - cI * TMP0;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV1P0_3( const cxtype F1[],   // input: wavefunction1[6]
                 const cxtype F2[],   // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2]);
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV2_0( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype V3[],   // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )     // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    const cxtype TMP1 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
                         F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    (*vertex) = COUP * - cI * TMP1;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV2_3( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
                         F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    const cxtype denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP2);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2] * OM3 * TMP2);
    V3[5] = denom * (-cI) * (-F1[2] * F2[4] - P3[3] * OM3 * TMP2 + F1[3] * F2[5]);
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV4_0( const cxtype F1[],  // input: wavefunction1[6]
               const cxtype F2[],  // input: wavefunction2[6]
               const cxtype V3[],  // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )    // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    const cxtype TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
                         F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    const cxtype TMP4 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
                         F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])));
    (*vertex) = COUP * (-one) * (+cI * (TMP3) + two * cI * (TMP4));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV4_3( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype TMP5 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
                         F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3])));
    const cxtype TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
                         F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    const cxtype denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-two * cI) * (OM3 * - half * P3[0] * (TMP2 + two * (TMP5)) +
                                  (+half * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * F2[3]));
    V3[3] = denom * (-two * cI) * (OM3 * - half * P3[1] * (TMP2 + two * (TMP5)) +
                                  (-half * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * F2[2]));
    V3[4] = denom * two * cI * (OM3 * half * P3[2] * (TMP2 + two * (TMP5)) +
                               (+half * cI * (F1[2] * F2[5]) - half * cI * (F1[3] * F2[4]) - cI *
                                (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
    V3[5] = denom * two * cI * (OM3 * half * P3[3] * (TMP2 + two * (TMP5)) +
                               (+half * (F1[2] * F2[4]) - half * (F1[3] * F2[5]) - F1[4] * F2[2] + F1[5] * F2[3]));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV2_4_0( const cxtype F1[],    // input: wavefunction1[6]
                 const cxtype F2[],    // input: wavefunction2[6]
                 const cxtype V3[],    // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype* vertex )      // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    // Note: inverting the order and computing TMP4 before TMP3 has increased C++ speed by ~1%
    const cxtype TMP4 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
                         F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])));
    const cxtype TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
                         F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    (*vertex) = ( -one ) * (COUP2 * (+cI * (TMP3) + two * cI * (TMP4)) + cI * (TMP3 * COUP1));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__
  void FFV2_4_3( const cxtype F1[],   // input: wavefunction1[6]
                 const cxtype F2[],   // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    const fptype P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    // Note: inverting the order and computing TMP2 before TMP5 has increased C++ speed by ~0.5%
    const cxtype TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
                         F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    const cxtype TMP5 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
                         F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3])));
    const cxtype denom = one/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-two * cI) *
      (COUP2 * (OM3 * - half * P3[0] * (TMP2 + two * (TMP5)) +
                (+half * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * F2[3])) +
       half * (COUP1 * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2)));
    V3[3] = denom * (-two * cI) *
      (COUP2 * (OM3 * - half * P3[1] * (TMP2 + two * (TMP5)) +
                (-half * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * F2[2])) -
       half * (COUP1 * (F1[2] * F2[5] + F1[3] * F2[4] + P3[1] * OM3 * TMP2)));
    V3[4] = denom * cI *
      (COUP2 * (OM3 * P3[2] * (TMP2 + two * (TMP5)) +
                (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) - two * cI * (F1[4] * F2[3]) + two * cI * (F1[5] * F2[2]))) +
       COUP1 * (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) + P3[2] * OM3 * TMP2));
    V3[5] = denom * two * cI *
      (COUP2 * (OM3 * half * P3[3] * (TMP2 + two * (TMP5)) +
                (+half * (F1[2] * F2[4]) - half * (F1[3] * F2[5]) - F1[4] * F2[2] + F1[5] * F2[3])) +
       half * (COUP1 * (F1[2] * F2[4] + P3[3] * OM3 * TMP2 - F1[3] * F2[5])));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

}  // end namespace
