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
  template<typename T_Acc>
  ALPAKA_FN_ACC
  void getMomentaInitial::operator()( T_Acc const &acc,
                                      const fptype energy, // input: energy
                                      fptype momenta1d[]   // output: momenta as AOSOA[npagM][npar][4][neppM]
                                    ) const
  {
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;
    {
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid
      const int ievt = idim;
      //printf( "getMomentaInitial: ievt %d\n", ievt );

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
  template<typename T_Acc>
  ALPAKA_FN_ACC
  void getMomentaFinal::operator()( T_Acc const &acc,
                                    const fptype energy,      // input: energy
                                    const fptype rnarray1d[], // input: random numbers in [0,1] as AOSOA[npagR][nparf][4][neppR]
                                    fptype momenta1d[],       // output: momenta as AOSOA[npagM][npar][4][neppM]
                                    fptype wgts[]             // output: weights[nevt]
                                  ) const
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

    const int neppR = mgOnGpu::neppR; // ASA layout: constant at compile-time
    fptype (*rnarray)[nparf][np4][neppR] = (fptype (*)[nparf][np4][neppR]) rnarray1d; // cast to multiD array pointer (AOSOA)
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    
    // initialization step: factorials for the phase space weight
    const fptype twopi = 8. * atan(acc, 1.);
    const fptype po2log = log(acc, twopi / 4.);
    fptype z[nparf];
    z[1] = po2log;
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = z[kpar - 1] + po2log - 2. * log(acc, fptype(kpar - 1));
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = (z[kpar] - log(acc, fptype(kpar)));

    {
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid
      const int ievt = idim;
      //printf( "getMomentaFinal:   ievt %d\n", ievt );

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
        const fptype s = sqrt(acc, 1. - c * c);
        const fptype f = twopi * r2;
        q[iparf][0] = -log(acc, r3 * r4);
        q[iparf][3] = q[iparf][0] * c;
        q[iparf][2] = q[iparf][0] * s * cos(acc, f);
        q[iparf][1] = q[iparf][0] * s * sin(acc, f);
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
      const fptype rmas = sqrt(acc, pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
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
        wt = (2. * nparf - 4.) * log(acc, energy) + z[nparf-1];

      // return for weighted massless momenta
      // nothing else to do in this event if all particles are massless (nm==0)

    }
    // ** END LOOP ON IEVT **

    return;
  }

  //--------------------------------------------------------------------------

  // Create and initialise a generator
  void createGenerator( mgGenerator_t* pgen )
  {
    *pgen = new mgGenerator;
#if defined MGONGPU_RANDTYPE_CURAND
    {
      // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
      const curandRngType_t type = CURAND_RNG_PSEUDO_MTGP32;          // 0.00082s | 0.00064s (FOR FAST TESTS)
      //const curandRngType_t type = CURAND_RNG_PSEUDO_XORWOW;        // 0.049s   | 0.0016s 
      //const curandRngType_t type = CURAND_RNG_PSEUDO_MRG32K3A;      // 0.71s    | 0.0012s  (better but slower, especially in c++)
      //const curandRngType_t type = CURAND_RNG_PSEUDO_MT19937;       // 21s      | 0.021s
      //const curandRngType_t type = CURAND_RNG_PSEUDO_PHILOX4_32_10; // 0.024s   | 0.00026s (used to segfault?)
#if defined MGONGPU_RAND_ONDEVICE
      checkRand( curandCreateGenerator( &(*pgen)->gen, type ) );
#elif defined MGONGPU_RAND_ONHOST
      checkRand( curandCreateGeneratorHost( &(*pgen)->gen, type ) );
#endif
      //checkRand( curandSetGeneratorOrdering( (*pgen)->gen, CURAND_ORDERING_PSEUDO_LEGACY ) ); // CUDA 11
      checkRand( curandSetGeneratorOrdering( (*pgen)->gen, CURAND_ORDERING_PSEUDO_BEST ) );
    }
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
    {
      const alsimple::alsrandType_t type = alsimple::ALSIMPLE_RNG_PSEUDO_MT19937;
#if defined MGONGPU_RAND_ONHOST
      checkRand( alsimple::createGeneratorHost( &(*pgen)->gen, type ) );
#elif defined MGONGPU_RAND_ONDEVICE
#error random number generation on device not currently supported with alsimple
#endif
    }
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif
  }

  //--------------------------------------------------------------------------

  // Seed a generator
  void seedGenerator( mgGenerator_t gen, unsigned long long seed )
  {
#if defined MGONGPU_RANDTYPE_CURAND
    checkRand( curandSetPseudoRandomGeneratorSeed( gen->gen, seed ) );
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
    checkRand( alsimple::randSetPseudoRandomGeneratorSeed( gen->gen, seed ) );
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif
  }

  //--------------------------------------------------------------------------

  // Destroy a curand generator
  void destroyGenerator( mgGenerator_t gen )
  {
#if defined MGONGPU_RANDTYPE_CURAND
    checkRand( curandDestroyGenerator( gen->gen ) );
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
    checkRand( alsimple::randDestroyGenerator( gen->gen ) );
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif
    delete gen;
  }

  //--------------------------------------------------------------------------

  // Bulk-generate the random numbers needed to process nevt events in rambo
  // ** NB: the random numbers are always produced in the same order and are interpreted as an AOSOA
  // AOSOA: rnarray[npagR][nparf][np4][neppR] where nevt=npagR*neppR
  void generateRnarray( mgGenerator_t gen,     // input: generator
                        fptype rnarray1d[],    // input: random numbers in [0,1] as AOSOA[npagR][nparf][4][neppR]
                        const int nevt )       // input: #events
  {
#if defined MGONGPU_RANDTYPE_CURAND
#if defined MGONGPU_FPTYPE_DOUBLE
    checkRand( curandGenerateUniformDouble( gen->gen, rnarray1d, np4*nparf*nevt ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkRand( curandGenerateUniform( gen->gen, rnarray1d, np4*nparf*nevt ) );
#endif
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
#if defined MGONGPU_FPTYPE_DOUBLE
    checkRand( alsimple::randGenerateUniformDouble( gen->gen, rnarray1d, np4*nparf*nevt ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkRand( alsimple::randGenerateUniform( gen->gen, rnarray1d, np4*nparf*nevt ) );
#endif
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif
  }

  //--------------------------------------------------------------------------

}
