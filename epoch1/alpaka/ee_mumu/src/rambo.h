#ifndef RAMBO_H
#define RAMBO_H 1

#include <cassert>
#include <cuda_to_cupla.hpp>

// using namespace alpaka;

#if defined MGONGPU_RANDTYPE_CURAND
#include "curand.h"
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
#include "alsimple.h"
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

//--------------------------------------------------------------------------

#define checkRand( code )                     \
  { assertRand( code, __FILE__, __LINE__ ); }

#if defined MGONGPU_RANDTYPE_CURAND
inline void assertRand( curandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != CURAND_STATUS_SUCCESS )
  {
    printf( "assertRand: %s %d\n", file, line );
    if ( abort ) assert( code == CURAND_STATUS_SUCCESS );
  }
}
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
inline void assertRand( alsimple::randStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != 0 )
  {
    printf( "assertRand: %s %d\n", file, line );
    if ( abort ) assert( code == 0 );
  }
}
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif

//--------------------------------------------------------------------------

class mgGenerator {
public:
#if defined MGONGPU_RANDTYPE_CURAND
  curandGenerator_t gen;
#elif defined MGONGPU_RANDTYPE_ALSIMPLE
  alsimple::randGenerator_t gen;
#else
#error unrecognised MGONGPU_RANDTYPE selection
#endif
};

typedef mgGenerator* mgGenerator_t;

//--------------------------------------------------------------------------

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace rambo2toNm0
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
    void operator()( T_Acc const &acc,
                     const fptype energy,       // input: energy
                     fptype momenta1d[]) const; // output: momenta as AOSOA[npagM][npar][4][neppM]
  };

  //--------------------------------------------------------------------------

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  struct getMomentaFinal {
    template<typename T_Acc>
    ALPAKA_FN_ACC
    void operator()( T_Acc const &acc,
                     const fptype energy,      // input: energy
                     const fptype rnarray1d[], // input: random numbers in [0,1] as AOSOA[npagR][nparf][4][neppR]
                     fptype momenta1d[],       // output: momenta as AOSOA[npagM][npar][4][neppM]
                     fptype wgts[]             // output: weights[nevt]
                    ) const;
  };

  //--------------------------------------------------------------------------

  // Create and initialise a generator
  void createGenerator( mgGenerator_t* pgen );

  //--------------------------------------------------------------------------
  
  // Seed a generator
  void seedGenerator( mgGenerator_t gen, unsigned long long seed );
  
  //--------------------------------------------------------------------------

  // Destroy a generator
  void destroyGenerator( mgGenerator_t gen );

  //--------------------------------------------------------------------------

  // Bulk-generate the random numbers needed to process nevt events in rambo
  // ** NB: the random numbers are always produced in the same order and are interpreted as an AOSOA
  // AOSOA: rnarray[npagR][nparf][np4][neppR] where nevt=npagR*neppR
  void generateRnarray( mgGenerator_t gen,     // input: curand generator
                        fptype rnarray1d[],    // input: random numbers in [0,1] as AOSOA[npagR][nparf][4][neppR]
                        const int nevt );      // input: #events

  //--------------------------------------------------------------------------

}

#include "rambo.cc"

#endif // RAMBO_H 1
