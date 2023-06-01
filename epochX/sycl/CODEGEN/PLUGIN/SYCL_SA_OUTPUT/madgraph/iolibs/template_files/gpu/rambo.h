// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
//==========================================================================

#ifndef RAMBO_H
#define RAMBO_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

//--------------------------------------------------------------------------


//--------------------------------------------------------------------------

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace rambo2toNm0
{

  static constexpr int np4 = mgOnGpu::np4;
  static constexpr int npari = mgOnGpu::npari;
  static constexpr int nparf = mgOnGpu::nparf;
  static constexpr int npar = mgOnGpu::npar;

  //--------------------------------------------------------------------------
  // Memory access functions
  SYCL_EXTERNAL
  inline fptype& kernelAccessIp4IparIevt(
          fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
          const int ip4,
          const int ipar,
          const size_t ievt
          ) {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    static constexpr int np4 = mgOnGpu::np4;
    static constexpr int npar = mgOnGpu::npar;
    static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    return momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  SYCL_EXTERNAL
  inline const fptype& kernelAccessIp4IparfIevt(
          const fptype* buffer, 
          const int ip4,
          const int iparf,
          const size_t ievt
          ) {
    static constexpr int neppR = 8; // FIXME HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
    const int ipagR = ievt/neppR; // #event "R-page"
    const int ieppR = ievt%neppR; // #event in the current event R-page
    return buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR]; // AOSOA[ipagR][iparf][ip4][ieppR]
  }

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  SYCL_EXTERNAL
  void ramboGetMomentaInitial( const fptype energy,  // input: energy
                               fptype* momenta,      // output: momenta for one event or for a set of events
                               const size_t ievt )
  {
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;
    kernelAccessIp4IparIevt( momenta, 0, 0, ievt ) = energy1;
    kernelAccessIp4IparIevt( momenta, 1, 0, ievt ) = 0;
    kernelAccessIp4IparIevt( momenta, 2, 0, ievt ) = 0;
    kernelAccessIp4IparIevt( momenta, 3, 0, ievt ) = mom;
    kernelAccessIp4IparIevt( momenta, 0, 1, ievt ) = energy2;
    kernelAccessIp4IparIevt( momenta, 1, 1, ievt ) = 0;
    kernelAccessIp4IparIevt( momenta, 2, 1, ievt ) = 0;
    kernelAccessIp4IparIevt( momenta, 3, 1, ievt ) = -mom;
  }

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  SYCL_EXTERNAL
  void ramboGetMomentaFinal( const fptype energy,      // input: energy
                             const fptype* rnarray,    // input: random numbers in [0,1] for one event or for a set of events
                             fptype* momenta,          // output: momenta for one event or for a set of events
                             fptype* wgts,             // output: weights for one event or for a set of events
                             const size_t ievt )
  {
    /****************************************************************************
     *                       rambo                                              *
     *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                       *
     *                                                                          *
     *    a democratic multi-particle phase space generator                     *
     *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                      *
     *    this is version 1.0 -  written by r. kleiss                           *
     *    -- adjusted by hans kuijf, weights are logarithmic (1990-08-20)       *
     *    -- adjusted by madgraph@sheffield_gpu_hackathon team (2020-07-29)     *
     *                                                                          *
     ****************************************************************************/

    // output weight
    fptype& wt = wgts[ievt];

    // AV special case nparf==1 (issue #358)
    if constexpr ( nparf == 1 )
    {
      bool first = true;
      if ( first )
      {
        first = false;
      }
      const int iparf = 0;
      for ( int i4 = 0; i4 < np4; i4++ )
      {
        kernelAccessIp4IparIevt( momenta, i4, iparf+npari, ievt ) = 0;
        for ( int ipari = 0; ipari < npari; ipari++ )
        {
          kernelAccessIp4IparIevt( momenta, i4, iparf+npari, ievt ) += kernelAccessIp4IparIevt( momenta, i4, ipari, ievt );
        }
      }
      wt = 1;
      return;
    }

    // initialization step: factorials for the phase space weight
    const fptype twopi = 8. * sycl::atan(1.);
    const fptype po2log = sycl::log(twopi / 4.);
    fptype z[nparf];
    z[1] = po2log;
    for ( int kpar = 2; kpar < nparf; kpar++ ) z[kpar] = z[kpar - 1] + po2log - 2. * sycl::log(fptype(kpar - 1));
    for ( int kpar = 2; kpar < nparf; kpar++ ) z[kpar] = (z[kpar] - sycl::log(fptype(kpar)));

    // generate n massless momenta in infinite phase space
    fptype q[nparf][np4];
    for ( int iparf = 0; iparf < nparf; iparf++ )
    {
      const fptype r1 = kernelAccessIp4IparfIevt( rnarray, 0, iparf, ievt );
      const fptype r2 = kernelAccessIp4IparfIevt( rnarray, 1, iparf, ievt );
      const fptype r3 = kernelAccessIp4IparfIevt( rnarray, 2, iparf, ievt );
      const fptype r4 = kernelAccessIp4IparfIevt( rnarray, 3, iparf, ievt );
      const fptype c = 2. * r1 - 1.;
      const fptype s = sycl::sqrt(1. - c * c);
      const fptype f = twopi * r2;
      q[iparf][0] = -sycl::log(r3 * r4);
      q[iparf][3] = q[iparf][0] * c;
      q[iparf][2] = q[iparf][0] * s * sycl::cos(f);
      q[iparf][1] = q[iparf][0] * s * sycl::sin(f);
    }

    // calculate the parameters of the conformal transformation
    fptype r[np4];
    fptype b[np4-1];
    for ( int i4 = 0; i4 < np4; i4++ ) r[i4] = 0.;
    for ( int iparf = 0; iparf < nparf; iparf++ )
    {
      for ( int i4 = 0; i4 < np4; i4++ ) r[i4] = r[i4] + q[iparf][i4];
    }
    const fptype rmas = sycl::sqrt(sycl::pown(r[0], 2) - sycl::pown(r[3], 2) - sycl::pown(r[2], 2) - sycl::pown(r[1], 2));
    for ( int i4 = 1; i4 < np4; i4++ ) b[i4-1] = -r[i4] / rmas;
    const fptype g = r[0] / rmas;
    const fptype a = 1. / (1. + g);
    const fptype x0 = energy / rmas;

    // transform the q's conformally into the p's (i.e. the 'momenta')
    for ( int iparf = 0; iparf < nparf; iparf++ )
    {
      fptype bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
      for ( int i4 = 1; i4 < np4; i4++ )
      {
        kernelAccessIp4IparIevt( momenta, i4, iparf+npari, ievt ) = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
      }
      kernelAccessIp4IparIevt( momenta, 0, iparf+npari, ievt ) = x0 * (g * q[iparf][0] + bq);
    }

    // calculate weight (NB return log of weight)
    wt = po2log;
    if ( nparf != 2 ) wt = (2. * nparf - 4.) * sycl::log(energy) + z[nparf-1];

    // return for weighted massless momenta
    // nothing else to do in this event if all particles are massless (nm==0)

    return;
  }

}

#endif // RAMBO_H 1
