#include <CL/sycl.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "rambo.h"

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace rambo2toNm0
{

  //--------------------------------------------------------------------------

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  SYCL_EXTERNAL void getMomentaInitial(
          const fptype energy,      // input: energy
          fptype momenta1d[],       // output: momenta as AOSOA[npagM][npar][4][neppM]
          sycl::nd_item<3> item_ct1 // SYCL item
          ) {
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;
    {
      const int idim = item_ct1.get_local_range().get(2) * item_ct1.get_group(2) + item_ct1.get_local_id(2); // event# == threadid
      const int ievt = idim;
      const int ipagM = ievt/neppM; // #eventpage in this iteration
      const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
      momenta[ipagM][0][0][ieppM] = energy1;
      momenta[ipagM][0][1][ieppM] = 0;
      momenta[ipagM][0][2][ieppM] = 0;
      momenta[ipagM][0][3][ieppM] = mom;
      momenta[ipagM][1][0][ieppM] = energy2;
      momenta[ipagM][1][1][ieppM] = 0;
      momenta[ipagM][1][2][ieppM] = 0;
      momenta[ipagM][1][3][ieppM] = -mom;
    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  SYCL_EXTERNAL void getMomentaFinal(
          const fptype energy,      // input: energy
          const fptype rnarray1d[], // input: random numbers in [0,1] as AOSOA[npagR][nparf][4][neppR]
          fptype momenta1d[],       // output: momenta as AOSOA[npagM][npar][4][neppM]
          fptype wgts[],            // output: weights[nevt]
          sycl::nd_item<3> item_ct1 // SYCL item
          ) {

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

    const int neppR = mgOnGpu::neppR; // ASA layout: constant at compile-time
    fptype (*rnarray)[nparf][np4][neppR] = (fptype (*)[nparf][np4][neppR]) rnarray1d; // cast to multiD array pointer (AOSOA)
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    
    // initialization step: factorials for the phase space weight
    const fptype twopi = 8. * sycl::atan(1.);
    const fptype po2log = sycl::log(twopi / 4.);
    fptype z[nparf];
    z[1] = po2log;
    for (int kpar = 2; kpar < nparf; kpar++)
        z[kpar] = z[kpar - 1] + po2log - 2. * sycl::log(fptype(kpar - 1));
    for (int kpar = 2; kpar < nparf; kpar++)
        z[kpar] = (z[kpar] - sycl::log(fptype(kpar)));

    {
        const int idim = item_ct1.get_local_range().get(2) * item_ct1.get_group(2) + item_ct1.get_local_id(2); // event# == threadid
        const int ievt = idim;

        const int ipagR = ievt/neppR; // #eventpage in this iteration
        const int ieppR = ievt%neppR; // #event in the current eventpage in this iteration
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration

        fptype& wt = wgts[ievt];

        // generate n massless momenta in infinite phase space
        fptype q[nparf][np4];
        for (int iparf = 0; iparf < nparf; iparf++) {
            const fptype r1 = rnarray[ipagR][iparf][0][ieppR];
            const fptype r2 = rnarray[ipagR][iparf][1][ieppR];
            const fptype r3 = rnarray[ipagR][iparf][2][ieppR];
            const fptype r4 = rnarray[ipagR][iparf][3][ieppR];
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
      for (int i4 = 0; i4 < np4; i4++)
          r[i4] = 0.;
      for (int iparf = 0; iparf < nparf; iparf++) {
          for (int i4 = 0; i4 < np4; i4++)
            r[i4] = r[i4] + q[iparf][i4];
      }
      const fptype rmas = sycl::sqrt(r[0] * r[0] - r[3] * r[3] - r[2] * r[2] - r[1] * r[1]);
      for (int i4 = 1; i4 < np4; i4++)
          b[i4-1] = -r[i4] / rmas;
      const fptype g = r[0] / rmas;
      const fptype a = 1. / (1. + g);
      const fptype x0 = energy / rmas;

      // transform the q's conformally into the p's (i.e. the 'momenta')
      for (int iparf = 0; iparf < nparf; iparf++) {
          fptype bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
          for (int i4 = 1; i4 < np4; i4++)
              momenta[ipagM][iparf+npari][i4][ieppM] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
          momenta[ipagM][iparf+npari][0][ieppM] = x0 * (g * q[iparf][0] + bq);
      }

      // calculate weight (NB return log of weight)
      wt = po2log;
      if (nparf != 2)
          wt = (2. * nparf - 4.) * sycl::log(energy) + z[nparf-1];

      // return for weighted massless momenta
      // nothing else to do in this event if all particles are massless (nm==0)

    } // ** END LOOP ON IEVT **
  }
}

