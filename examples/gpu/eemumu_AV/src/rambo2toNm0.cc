#include <cmath>
#include <cstdlib>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "Random.h"

#ifndef __CUDACC__
#include "rambo2toNm0.h"
#endif

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
#ifdef __CUDACC__
namespace grambo2toNm0
#else
namespace rambo2toNm0
#endif
{
  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)

  const int npari = 2; // #particles in the initial state
  const int nparf = 2; // #particles in the final state
  const int npar = npari + nparf; // #particles in total

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
#ifdef __CUDACC__
__global__
#endif
  void getMomentaInitial( const double energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                          double momenta1d[],     // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                          double momenta1d[],     // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                          double momenta1d[],     // output: momenta as AOS[nevt][npar][4]
#endif
                          const int nevt )        // input: #events
  {
#if defined MGONGPU_LAYOUT_ASA
    using mgOnGpu::nepp;
    double (*momenta)[npar][np4][nepp] = (double (*)[npar][np4][nepp]) momenta1d; // cast to multiD array pointer (AOSOA)
#elif defined MGONGPU_LAYOUT_SOA
    double (*momenta)[np4][nevt] = (double (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#elif defined MGONGPU_LAYOUT_AOS
    double (*momenta)[npar][np4] = (double (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif
    const double energy1 = energy/2;
    const double energy2 = energy/2;
    const double mom = energy/2;
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt) 
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // 0 to ndim-1
      const int ievt = idim;
      //printf( "getMomentaInitial: ievt %d\n", ievt );
#endif    
#if defined MGONGPU_LAYOUT_ASA
      const int ipag = ievt/nepp; // #eventpage in this iteration
      const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
      momenta[ipag][0][0][iepp] = energy1;
      momenta[ipag][0][1][iepp] = 0;
      momenta[ipag][0][2][iepp] = 0;
      momenta[ipag][0][3][iepp] = mom;
      momenta[ipag][1][0][iepp] = energy2;
      momenta[ipag][1][1][iepp] = 0;
      momenta[ipag][1][2][iepp] = 0;
      momenta[ipag][1][3][iepp] = -mom;
#elif defined MGONGPU_LAYOUT_SOA
      momenta[0][0][ievt] = energy1;
      momenta[0][1][ievt] = 0;
      momenta[0][2][ievt] = 0;
      momenta[0][3][ievt] = mom;
      momenta[1][0][ievt] = energy2;
      momenta[1][1][ievt] = 0;
      momenta[1][2][ievt] = 0;
      momenta[1][3][ievt] = -mom;
#elif defined MGONGPU_LAYOUT_AOS
      momenta[ievt][0][0] = energy1;
      momenta[ievt][0][1] = 0;
      momenta[ievt][0][2] = 0;
      momenta[ievt][0][3] = mom;
      momenta[ievt][1][0] = energy2;
      momenta[ievt][1][1] = 0;
      momenta[ievt][1][2] = 0;
      momenta[ievt][1][3] = -mom;
#endif
    }
    // ** END LOOP ON IEVT **
  }

#ifdef __CUDACC__
__global__
#endif
  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  void getMomentaFinal( const double energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                        const double rnarray1d[], // input: randomnumbers in [0,1] as AOSOA[npag][nparf][4][nepp]
                        double momenta1d[],       // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                        const double rnarray1d[], // input: randomnumbers in [0,1] as SOA[nparf][4][nevt]
                        double momenta1d[],       // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                        const double rnarray1d[], // input: randomnumbers in [0,1] as SOA[nevt][nparf][4]
                        double momenta1d[],       // output: momenta as AOS[nevt][npar][4]
#endif
                        double wgts[],            // output: weights[nevt]
                        const int nevt )          // input: #events
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
    const double twopi = 8. * atan(1.);
    const double po2log = log(twopi / 4.);
    double z[nparf];
    z[1] = po2log;
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = z[kpar - 1] + po2log - 2. * log(double(kpar - 1));
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = (z[kpar] - log(double(kpar)));

#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt) 
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // 0 to ndim-1
      const int ievt = idim;
      //printf( "getMomentaFinal:   ievt %d\n", ievt );
#endif    

#if defined MGONGPU_LAYOUT_ASA
      using mgOnGpu::nepp;
      const int ipag = ievt/nepp; // #eventpage in this iteration
      const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
      double (*rnarray)[nparf][np4][nepp] = (double (*)[nparf][np4][nepp]) rnarray1d; // cast to multiD array pointer (AOSOA)
      double (*momenta)[npar][np4][nepp] = (double (*)[npar][np4][nepp]) momenta1d; // cast to multiD array pointer (AOSOA)
#elif defined MGONGPU_LAYOUT_SOA
      double (*rnarray)[np4][nevt] = (double (*)[np4][nevt]) rnarray1d; // cast to multiD array pointer (SOA)
      double (*momenta)[np4][nevt] = (double (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#elif defined MGONGPU_LAYOUT_AOS
      double (*rnarray)[nparf][np4] = (double (*)[nparf][np4]) rnarray1d; // cast to multiD array pointer (AOS)
      double (*momenta)[npar][np4] = (double (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif
      double& wt = wgts[ievt];

      // generate n massless momenta in infinite phase space
      double q[nparf][np4];
      for (int iparf = 0; iparf < nparf; iparf++) {
#if defined MGONGPU_LAYOUT_ASA
        const double r1 = rnarray[ipag][iparf][0][iepp];
        const double r2 = rnarray[ipag][iparf][1][iepp];
        const double r3 = rnarray[ipag][iparf][2][iepp];
        const double r4 = rnarray[ipag][iparf][3][iepp];
#elif defined MGONGPU_LAYOUT_SOA
        const double r1 = rnarray[iparf][0][ievt];
        const double r2 = rnarray[iparf][1][ievt];
        const double r3 = rnarray[iparf][2][ievt];
        const double r4 = rnarray[iparf][3][ievt];
#elif defined MGONGPU_LAYOUT_AOS
        const double r1 = rnarray[ievt][iparf][0];
        const double r2 = rnarray[ievt][iparf][1];
        const double r3 = rnarray[ievt][iparf][2];
        const double r4 = rnarray[ievt][iparf][3];
#endif
        const double c = 2. * r1 - 1.;
        const double s = sqrt(1. - c * c);
        const double f = twopi * r2;
        q[iparf][0] = -log(r3 * r4);
        q[iparf][3] = q[iparf][0] * c;
        q[iparf][2] = q[iparf][0] * s * cos(f);
        q[iparf][1] = q[iparf][0] * s * sin(f);
      }

      // calculate the parameters of the conformal transformation
      double r[np4];
      double b[np4-1];
      for (int i4 = 0; i4 < np4; i4++)
        r[i4] = 0.;
      for (int iparf = 0; iparf < nparf; iparf++) {
        for (int i4 = 0; i4 < np4; i4++)
          r[i4] = r[i4] + q[iparf][i4];
      }
      const double rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
      for (int i4 = 1; i4 < np4; i4++)
        b[i4-1] = -r[i4] / rmas;
      const double g = r[0] / rmas;
      const double a = 1. / (1. + g);
      const double x0 = energy / rmas;

      // transform the q's conformally into the p's (i.e. the 'momenta')
      for (int iparf = 0; iparf < nparf; iparf++) {
        double bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
#if defined MGONGPU_LAYOUT_ASA
        for (int i4 = 1; i4 < np4; i4++)
          momenta[ipag][iparf+npari][i4][iepp] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta[ipag][iparf+npari][0][iepp] = x0 * (g * q[iparf][0] + bq);
#elif defined MGONGPU_LAYOUT_SOA
        for (int i4 = 1; i4 < np4; i4++)
          momenta[iparf+npari][i4][ievt] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta[iparf+npari][0][ievt] = x0 * (g * q[iparf][0] + bq);
#elif defined MGONGPU_LAYOUT_AOS
        for (int i4 = 1; i4 < np4; i4++)
          momenta[ievt][iparf+npari][i4] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta[ievt][iparf+npari][0] = x0 * (g * q[iparf][0] + bq);
#endif
      }

      // calculate weight (NB return log of weight)
      wt = po2log;
      if (nparf != 2)
        wt = (2. * nparf - 4.) * log(energy) + z[nparf-1];

#ifndef __CUDACC__
      // issue warnings if weight is too small or too large
      static int iwarn[5] = {0,0,0,0,0};
      if (wt < -180.) {
        if (iwarn[0] <= 5)
          std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
        iwarn[0] = iwarn[0] + 1;
      }
      if (wt > 174.) {
        if (iwarn[1] <= 5)
          std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
        iwarn[1] = iwarn[1] + 1;
      }
#endif

      // return for weighted massless momenta
      // nothing else to do in this event if all particles are massless (nm==0)

    }
    // ** END LOOP ON IEVT **

    return;
  }

  // Generate the random numbers needed to process nevt events in rambo
#if defined MGONGPU_LAYOUT_ASA
  // AOSOA: rnarray[npag][nparf][np4][nepp] where nevt=npag*nepp
#elif defined MGONGPU_LAYOUT_SOA
  // SOA: rnarray[nparf][np4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
  // AOS: rnarray[nevt][nparf][np4]
#endif
  void generateRnArray( double rnarray1d[], // output: randomnumbers in [0,1]
                        const int nevt,     // input: #events
                        const int iiter )   // input: iteration#
  {
    const int np4 = 4; // 4 random numbers (like the dimension of 4-momenta) are needed for each particle
#if defined MGONGPU_LAYOUT_ASA
    using mgOnGpu::nepp;
    double (*rnarray)[nparf][np4][nepp] = (double (*)[nparf][np4][nepp]) rnarray1d; // cast to multiD array pointer (AOSOA)
#elif defined MGONGPU_LAYOUT_SOA
    double (*rnarray)[np4][nevt] = (double (*)[np4][nevt]) rnarray1d; // cast to multiD array pointer (SOA)
#elif defined MGONGPU_LAYOUT_AOS
    double (*rnarray)[nparf][np4] = (double (*)[nparf][np4]) rnarray1d; // cast to multiD array pointer (AOS)
#endif
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt) 
    {
      Random random; // initialise with rmarin(ij,kl) where 0<=ij<=31328 and 0<=kl<=30081
      int ievt_ij = (ievt + nevt*iiter) % 31329;
      int ievt_kl = (ievt + nevt*iiter) % 30082;
      random.rmarin( ievt_ij, ievt_kl );
#if defined MGONGPU_LAYOUT_ASA
      const int ipag = ievt/nepp; // #eventpage in this iteration
      const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
      for (int iparf = 0; iparf < nparf; iparf++)
        for (int ip4 = 0; ip4 < np4; ++ip4)
          rnarray[ipag][iparf][ip4][iepp] = random.ranmar2(); // AOSOA[npag][nparf][np4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
      for (int iparf = 0; iparf < nparf; iparf++)
        for (int ip4 = 0; ip4 < np4; ++ip4)
          rnarray[iparf][ip4][ievt] = random.ranmar2(); // SOA[nparf][np4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
      for (int iparf = 0; iparf < nparf; iparf++)
        for (int ip4 = 0; ip4 < np4; ++ip4)
          rnarray[ievt][iparf][ip4] = random.ranmar2(); // AOS[nevt][nparf][np4]
#endif
    }
    // ** END LOOP ON IEVT **
    return;
  }

}
