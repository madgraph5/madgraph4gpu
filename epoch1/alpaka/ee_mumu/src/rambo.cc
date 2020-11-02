#include <cmath>
#include <cstdlib>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "rambo.h"

#include <iostream>

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace grambo2toNm0
{

  //--------------------------------------------------------------------------

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  template< typename T_Acc >
  ALPAKA_FN_ACC
  void getMomentaInitial::operator()( T_Acc const &acc,
                          const fptype energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                          fptype momenta1d[],     // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                          fptype momenta1d[],     // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                          fptype momenta1d[],     // output: momenta as AOS[nevt][npar][4]
#endif
                          const int nevt ) const  // input: #events
  {
#if defined MGONGPU_LAYOUT_ASA
    using mgOnGpu::nepp;
    fptype (*momenta)[npar][np4][nepp] = (fptype (*)[npar][np4][nepp]) momenta1d; // cast to multiD array pointer (AOSOA)
#elif defined MGONGPU_LAYOUT_SOA
    // Cast is impossible in CUDA C ("error: expression must have a constant value")
    //fptype (*momenta)[np4][nevt] = (fptype (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#elif defined MGONGPU_LAYOUT_AOS
    fptype (*momenta)[npar][np4] = (fptype (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;
    {
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // 0 to ndim-1
      const int ievt = idim;
      //printf( "getMomentaInitial: ievt %d\n", ievt );

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
      momenta1d[0*np4*nevt + 0*nevt + ievt] = energy1;
      momenta1d[0*np4*nevt + 1*nevt + ievt] = 0;
      momenta1d[0*np4*nevt + 2*nevt + ievt] = 0;
      momenta1d[0*np4*nevt + 3*nevt + ievt] = mom;
      momenta1d[1*np4*nevt + 0*nevt + ievt] = energy2;
      momenta1d[1*np4*nevt + 1*nevt + ievt] = 0;
      momenta1d[1*np4*nevt + 2*nevt + ievt] = 0;
      momenta1d[1*np4*nevt + 3*nevt + ievt] = -mom;
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

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  template< typename T_Acc >
  ALPAKA_FN_ACC
  void getMomentaFinal::operator()( T_Acc const &acc,
                        const fptype energy,      // input: energy
                        const fptype rnarray1d[], // input: randomnumbers in [0,1] as AOSOA[npag][nparf][4][nepp]
#if defined MGONGPU_LAYOUT_ASA
                        fptype momenta1d[],       // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                        fptype momenta1d[],       // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                        fptype momenta1d[],       // output: momenta as AOS[nevt][npar][4]
#endif
                        fptype wgts[],            // output: weights[nevt]
                        const int nevt ) const    // input: #events
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
    const fptype twopi = 8. * atan(acc, 1.);
    const fptype po2log = log(acc, twopi / 4.);
    fptype z[nparf];
    z[1] = po2log;
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = z[kpar - 1] + po2log - 2. * log(acc,fptype(kpar - 1));
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = (z[kpar] - log(acc,fptype(kpar)));

    {
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // 0 to ndim-1
      const int ievt = idim;
      //printf( "getMomentaFinal:   ievt %d\n", ievt );

      using mgOnGpu::nepp;
      const int ipag = ievt/nepp; // #eventpage in this iteration
      const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
//dhs
//      fptype (*rnarray)[nparf][np4][nepp] = (fptype (*)[nparf][np4][nepp]) rnarray1d; // cast to multiD array pointer (AOSOA)
      fptype (*rnarray)[nepp][nparf][np4] = (fptype (*)[nepp][nparf][np4]) rnarray1d; // cast to multiD array pointer (debug aos for consistency with standalone)
#if defined MGONGPU_LAYOUT_ASA
      fptype (*momenta)[npar][np4][nepp] = (fptype (*)[npar][np4][nepp]) momenta1d; // cast to multiD array pointer (AOSOA)
#elif defined MGONGPU_LAYOUT_SOA
      // Cast is impossible in CUDA C ("error: expression must have a constant value")
      //fptype (*momenta)[np4][nevt] = (fptype (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#elif defined MGONGPU_LAYOUT_AOS
      fptype (*momenta)[npar][np4] = (fptype (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif
      fptype& wt = wgts[ievt];

      // generate n massless momenta in infinite phase space
      fptype q[nparf][np4];
      for (int iparf = 0; iparf < nparf; iparf++) {
//dhs
        const fptype r1 = rnarray[ipag][iepp][iparf][0];
        const fptype r2 = rnarray[ipag][iepp][iparf][1];
        const fptype r3 = rnarray[ipag][iepp][iparf][2];
        const fptype r4 = rnarray[ipag][iepp][iparf][3];

//        const fptype r1 = rnarray[ipag][iparf][0][iepp];
//        const fptype r2 = rnarray[ipag][iparf][1][iepp];
//        const fptype r3 = rnarray[ipag][iparf][2][iepp];
//        const fptype r4 = rnarray[ipag][iparf][3][iepp];
//dhs
//#include <mutex>
//static std::mutex mmtex;
//std::lock_guard<std::mutex> mgu(mmtex);
//std::cout << "iparf=" << iparf << " r1=" << r1 << " r2=" << r2 << " r3=" << r3 << " r4=" << r4 << std::endl;
        const fptype c = 2. * r1 - 1.;
        const fptype s = sqrt(acc, 1. - c * c);
        const fptype f = twopi * r2;
        q[iparf][0] = -log(acc,r3 * r4);
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
#if defined MGONGPU_LAYOUT_ASA
        for (int i4 = 1; i4 < np4; i4++)
          momenta[ipag][iparf+npari][i4][iepp] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta[ipag][iparf+npari][0][iepp] = x0 * (g * q[iparf][0] + bq);
#elif defined MGONGPU_LAYOUT_SOA
        for (int i4 = 1; i4 < np4; i4++)
          momenta1d[(iparf+npari)*np4*nevt + i4*nevt + ievt] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta1d[(iparf+npari)*np4*nevt + 0*nevt + ievt] = x0 * (g * q[iparf][0] + bq);
#elif defined MGONGPU_LAYOUT_AOS
        for (int i4 = 1; i4 < np4; i4++)
          momenta[ievt][iparf+npari][i4] = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        momenta[ievt][iparf+npari][0] = x0 * (g * q[iparf][0] + bq);
#endif
      }

      // calculate weight (NB return log of weight)
      wt = po2log;
      if (nparf != 2)
        wt = (2. * nparf - 4.) * log(acc,energy) + z[nparf-1];

      // return for weighted massless momenta
      // nothing else to do in this event if all particles are massless (nm==0)

    }
    // ** END LOOP ON IEVT **

    return;
  }

  //--------------------------------------------------------------------------

  // Create and initialise a curand generator
  void createGenerator( mgRandGenerator_t* pgen )
  {
#if defined MGONGPU_CURAND_ONDEVICE || defined MGONGPU_CURAND_ONHOST
    // [NB Timings are for host generation of 32*256*1 events: rn(0) is 0.0012s]
    const curandRngType_t type = CURAND_RNG_PSEUDO_MTGP32;          // 0.0021s (FOR FAST TESTS)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_XORWOW;        // 1.13s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MRG32K3A;      // 10.5s (better but slower)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MT19937;       // 43s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_PHILOX4_32_10; // segfaults
#endif

#if defined MGONGPU_CURAND_ONDEVICE
    checkCurand( curandCreateGenerator( &pgen->rnGen, type ) );
#elif defined MGONGPU_CURAND_ONHOST
    checkCurand( curandCreateGeneratorHost( &pgen->rnGen, type ) );
#elif defined MGONGPU_RANDEXTRAS_ONHOST
    extras::extrasGenerator_t::rnCreate(&pgen->rnGen);
#endif

#if defined MGONGPU_CURAND_ONDEVICE || defined MGONGPU_CURAND_ONHOST
    //checkCurand( curandSetGeneratorOrdering( pgen->rnGen, CURAND_ORDERING_PSEUDO_LEGACY ) ); // CUDA 11
    checkCurand( curandSetGeneratorOrdering( pgen->rnGen, CURAND_ORDERING_PSEUDO_BEST ) );
#endif
  }

  //--------------------------------------------------------------------------

  // Seed a curand generator
  void seedGenerator( mgRandGenerator_t gen, unsigned long long seed )
  {
#if defined MGONGPU_RANDEXTRAS_ONHOST
    gen.rnGen->setseed(seed);
#elif defined MGONGPU_CURAND_ONDEVICE || defined MGONGPU_CURAND_ONHOST
    checkCurand( curandSetPseudoRandomGeneratorSeed( gen.rnGen, seed ) );
#endif
  }

  //--------------------------------------------------------------------------

  // Destroy a curand generator
  void destroyGenerator( mgRandGenerator_t gen )
  {
#if defined MGONGPU_RANDEXTRAS_ONHOST
    extras::extrasGenerator_t::rnDestroy(gen.rnGen);
#elif defined MGONGPU_CURAND_ONDEVICE || defined MGONGPU_CURAND_ONHOST
    checkCurand( curandDestroyGenerator( gen.rnGen ) );
#endif
  }

  //--------------------------------------------------------------------------

  // Bulk-generate (using curand) the random numbers needed to process nevt events in rambo
  // ** NB: the random numbers are always produced in the same order and are interpreted as an AOSOA
  // AOSOA: rnarray[npag][nparf][np4][nepp] where nevt=npag*nepp
  void generateRnArray( mgRandGenerator_t gen, // input: curand generator
                        fptype rnarray1d[],    // output: randomnumbers in [0,1]
                        const int nevt )       // input: #events
  {
#if defined MGONGPU_RANDEXTRAS_ONHOST
    gen.rnGen->uniform(rnarray1d, np4*nparf*nevt);
#else
#if defined MGONGPU_FPTYPE_DOUBLE
    checkCurand( curandGenerateUniformDouble( gen.rnGen, rnarray1d, np4*nparf*nevt ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkCurand( curandGenerateUniform( gen.rnGen, rnarray1d, np4*nparf*nevt ) );
#endif
#endif
  }

  //--------------------------------------------------------------------------

}
