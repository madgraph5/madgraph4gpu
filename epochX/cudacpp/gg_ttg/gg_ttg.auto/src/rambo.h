#include "mgOnGpuConfig.h"
#include "mgOnGpuFptypes.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  using mgOnGpu::np4;
  using mgOnGpu::npari;
  using mgOnGpu::nparf;
  using mgOnGpu::npar;

  //--------------------------------------------------------------------------

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  template<class M_ACCESS>
  __host__ __device__
  void ramboGetMomentaInitial( const fptype energy,  // input: energy
                               fptype* momenta )     // output: momenta for one event or for a set of events
  {
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 0, 0 ) = energy1;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 1, 0 ) = 0;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 2, 0 ) = 0;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 3, 0 ) = mom;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 0, 1 ) = energy2;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 1, 1 ) = 0;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 2, 1 ) = 0;
    M_ACCESS::kernelAccessIp4Ipar( momenta, 3, 1 ) = -mom;
  }

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  template<class R_ACCESS, class M_ACCESS, class W_ACCESS>
  __host__ __device__
  void ramboGetMomentaFinal( const fptype energy,      // input: energy
                             const fptype* rnarray,    // input: random numbers in [0,1] for one event or for a set of events
                             fptype* momenta,          // output: momenta for one event or for a set of events
                             fptype* wgts )            // output: weights for one event or for a set of events
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

    // initialization step: factorials for the phase space weight
    const fptype twopi = 8. * atan(1.);
    const fptype po2log = log(twopi / 4.);
    fptype z[nparf];
    z[1] = po2log;
    for ( int kpar = 2; kpar < nparf; kpar++ ) z[kpar] = z[kpar - 1] + po2log - 2. * log(fptype(kpar - 1));
    for ( int kpar = 2; kpar < nparf; kpar++ ) z[kpar] = (z[kpar] - log(fptype(kpar)));

    // output weight
    fptype& wt = W_ACCESS::kernelAccess( wgts );

    // generate n massless momenta in infinite phase space
    fptype q[nparf][np4];
    for ( int iparf = 0; iparf < nparf; iparf++ )
    {
      const fptype r1 = R_ACCESS::kernelAccessIp4IparfConst( rnarray, 0, iparf );
      const fptype r2 = R_ACCESS::kernelAccessIp4IparfConst( rnarray, 1, iparf );
      const fptype r3 = R_ACCESS::kernelAccessIp4IparfConst( rnarray, 2, iparf );
      const fptype r4 = R_ACCESS::kernelAccessIp4IparfConst( rnarray, 3, iparf );
      const fptype c = 2. * r1 - 1.;
      const fptype s = sqrt(1. - c * c);
      const fptype f = twopi * r2;
      q[iparf][0] = -log(r3 * r4);
      q[iparf][3] = q[iparf][0] * c;
      q[iparf][2] = q[iparf][0] * s * cos(f);
      q[iparf][1] = q[iparf][0] * s * sin(f);
    }

    // calculate the parameters of the conformal transformation
    fptype r[np4];
    fptype b[np4-1];
    for ( int i4 = 0; i4 < np4; i4++ ) r[i4] = 0.;
    for ( int iparf = 0; iparf < nparf; iparf++ )
    {
      for ( int i4 = 0; i4 < np4; i4++ ) r[i4] = r[i4] + q[iparf][i4];
    }
    const fptype rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
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
        M_ACCESS::kernelAccessIp4Ipar( momenta, i4, iparf+npari ) = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
      }
      M_ACCESS::kernelAccessIp4Ipar( momenta, 0, iparf+npari ) = x0 * (g * q[iparf][0] + bq);
    }

    // calculate weight (NB return log of weight)
    wt = po2log;
    if ( nparf != 2 ) wt = (2. * nparf - 4.) * log(energy) + z[nparf-1];

#ifndef __CUDACC__
    // issue warnings if weight is too small or too large
    static int iwarn[5] = {0,0,0,0,0};
    if ( wt < -180. )
    {
      if ( iwarn[0] <= 5 ) std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
      iwarn[0] = iwarn[0] + 1;
    }
    if ( wt > 174. )
    {
      if ( iwarn[1] <= 5 ) std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
      iwarn[1] = iwarn[1] + 1;
    }
#endif

    // return for weighted massless momenta
    // nothing else to do in this event if all particles are massless (nm==0)

    return;
  }

  //--------------------------------------------------------------------------

}

