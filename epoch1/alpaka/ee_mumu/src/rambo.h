#ifndef RAMBO2TONM0_H
#define RAMBO2TONM0_H 1

#include <cuda_to_cupla.hpp>
#include <cassert>

using namespace alpaka;

#if defined MGONGPU_RANDEXTRAS_ONHOST
#include "extras.h"
#elif defined MGONGPU_CURAND_ONHOST || defined MGONGPU_CURAND_ONDEVICE
#include "curand.h"
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

class mgRandGenerator_t {
public:
#if defined MGONGPU_RANDEXTRAS_ONHOST
  extras::extrasGenerator_t *rnGen;
#else
  curandGenerator_t rnGen;
#endif
};

//--------------------------------------------------------------------------

#if defined MGONGPU_CURAND_ONHOST || defined MGONGPU_CURAND_ONDEVICE
#define checkCurand( code )                     \
  { assertCurand( code, __FILE__, __LINE__ ); }

inline void assertCurand( curandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != CURAND_STATUS_SUCCESS )
  {
    printf( "CurandAssert: %s %d\n", file, line );
    if ( abort ) assert( code == CURAND_STATUS_SUCCESS );
  }
}
#endif

//--------------------------------------------------------------------------

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace grambo2toNm0
{

  using mgOnGpu::np4;
  using mgOnGpu::npari;
  using mgOnGpu::nparf;
  using mgOnGpu::npar;

  //--------------------------------------------------------------------------

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  struct getMomentaInitial {
    template<typename T_Acc>
    ALPAKA_FN_ACC
    void operator()(T_Acc const &acc,
                          const fptype energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                          fptype momenta1d[],     // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                          fptype momenta1d[],     // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                          fptype momenta1d[],     // output: momenta as AOS[nevt][npar][4]
#endif
                          const int nevt ) const; // input: #events
  };

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  struct getMomentaFinal {
    template<typename T_Acc>
    ALPAKA_FN_ACC
    void operator()(T_Acc const &acc,
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
                        const int nevt ) const;   // input: #events
  };

  //--------------------------------------------------------------------------

  // Create and initialise a curand generator
  void createGenerator( mgRandGenerator_t* pgen );

  //--------------------------------------------------------------------------
  
  // Seed a curand generator
  void seedGenerator( mgRandGenerator_t gen, unsigned long long seed );
  
  //--------------------------------------------------------------------------

  // Destroy a curand generator
  void destroyGenerator( mgRandGenerator_t gen );

  //--------------------------------------------------------------------------

  // Bulk-generate (using curand) the random numbers needed to process nevt events in rambo
  // ** NB: the random numbers are always produced in the same order and are interpreted as an AOSOA
  // AOSOA: rnarray[npag][nparf][np4][nepp] where nevt=npag*nepp
  void generateRnArray( mgRandGenerator_t gen, // input: curand generator
                        fptype rnarray1d[],    // output: randomnumbers in [0,1]
                        const int nevt );      // input: #events

  //--------------------------------------------------------------------------

}

#endif // RAMBO2TONM0_H 1
