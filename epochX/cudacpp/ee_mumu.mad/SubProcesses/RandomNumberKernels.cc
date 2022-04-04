#include "RandomNumberKernels.h"

#include "CommonRandomNumbers.h"
#include "CudaRuntime.h"
#include "MemoryBuffers.h"

#include <cassert>

#ifndef MGONGPU_HAS_NO_CURAND /* clang-format off */
#define checkCurand( code ){ assertCurand( code, __FILE__, __LINE__ ); }
inline void assertCurand( curandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != CURAND_STATUS_SUCCESS )
  {
    printf( "CurandAssert: %s %d\n", file, line );
    if ( abort ) assert( code == CURAND_STATUS_SUCCESS );
  }
}
#endif /* clang-format on */

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  CommonRandomNumberKernel::CommonRandomNumberKernel( BufferRandomNumbers& rnarray )
    : RandomNumberKernelBase( rnarray )
    , m_seed( 20211220 )
  {
    if( m_rnarray.isOnDevice() )
      throw std::runtime_error( "CommonRandomNumberKernel on host with a device random number array" );
  }

  //--------------------------------------------------------------------------

  void CommonRandomNumberKernel::generateRnarray()
  {
    std::vector<double> rnd = CommonRandomNumbers::generate<double>( m_rnarray.size(), m_seed ); // NB: generate as double (HARDCODED)
    std::copy( rnd.begin(), rnd.end(), m_rnarray.data() );                                       // NB: copy may imply a double-to-float conversion
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HAS_NO_CURAND
  CurandRandomNumberKernel::CurandRandomNumberKernel( BufferRandomNumbers& rnarray, const bool onDevice )
    : RandomNumberKernelBase( rnarray )
    , m_isOnDevice( onDevice )
  {
    if( m_isOnDevice )
    {
#ifdef __CUDACC__
      if( !m_rnarray.isOnDevice() )
        throw std::runtime_error( "CurandRandomNumberKernel on device with a host random number array" );
#else
      throw std::runtime_error( "CurandRandomNumberKernel does not support CurandDevice on CPU host" );
#endif
    }
    else
    {
      if( m_rnarray.isOnDevice() )
        throw std::runtime_error( "CurandRandomNumberKernel on host with a device random number array" );
    }
    // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
    const curandRngType_t type = CURAND_RNG_PSEUDO_MTGP32; //          0.00082s | 0.00064s (FOR FAST TESTS)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_XORWOW;        // 0.049s   | 0.0016s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MRG32K3A;      // 0.71s    | 0.0012s  (better but slower, especially in c++)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MT19937;       // 21s      | 0.021s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_PHILOX4_32_10; // 0.024s   | 0.00026s (used to segfault?)
    if( m_isOnDevice )
    {
      checkCurand( curandCreateGenerator( &m_rnGen, type ) );
    }
    else
    {
      checkCurand( curandCreateGeneratorHost( &m_rnGen, type ) );
    }
    //checkCurand( curandSetGeneratorOrdering( *&m_rnGen, CURAND_ORDERING_PSEUDO_LEGACY ) ); // CUDA 11
    checkCurand( curandSetGeneratorOrdering( *&m_rnGen, CURAND_ORDERING_PSEUDO_BEST ) );
  }

  //--------------------------------------------------------------------------

  CurandRandomNumberKernel::~CurandRandomNumberKernel()
  {
    checkCurand( curandDestroyGenerator( m_rnGen ) );
  }

  //--------------------------------------------------------------------------

  void CurandRandomNumberKernel::seedGenerator( const int seed )
  {
    //printf( "seedGenerator: seed %lld\n", seed );
    checkCurand( curandSetPseudoRandomGeneratorSeed( m_rnGen, seed ) );
  }

  //--------------------------------------------------------------------------

  void CurandRandomNumberKernel::generateRnarray()
  {
#if defined MGONGPU_FPTYPE_DOUBLE
    checkCurand( curandGenerateUniformDouble( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkCurand( curandGenerateUniform( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#endif
    //for ( int i=0; i<8; i++ ) printf("%f %f %f %f\n",m_rnarray.data()[i*4],m_rnarray.data()[i*4+2],m_rnarray.data()[i*4+2],m_rnarray.data()[i*4+3]);
  }

  //--------------------------------------------------------------------------
#endif
}
