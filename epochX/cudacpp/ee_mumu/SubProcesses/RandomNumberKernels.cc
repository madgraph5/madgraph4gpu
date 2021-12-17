#include "RandomNumberKernels.h"

#include "checkCuda.h"
#include "CommonRandomNumbers.h"

#include <cassert>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  RandomNumberKernelBase::RandomNumberKernelBase( const int nevt, const bool useHstRnarray )
    : m_nevt( nevt )
    , m_seed( -1 )
    , m_hstRnarray( nullptr )
    , m_devRnarray( nullptr )
  {
    // Allocate the host buffer for the output random numbers
    if ( useHstRnarray ) checkCuda( cudaMallocHost( &m_hstRnarray, nRnarray() * sizeof(fptype) ) );
    // Allocate the device buffer for the output random numbers
    checkCuda( cudaMalloc( &m_devRnarray, nRnarray() * sizeof(fptype) ) );
  }
#else
  RandomNumberKernelBase::RandomNumberKernelBase( const int nevt )
    : m_nevt( nevt )
    , m_seed( -1 )
    , m_hstRnarray( nullptr )
  {
    // Allocate the host buffer for the output random numbers
    m_hstRnarray = new( std::align_val_t{ cppAlign } ) fptype[ nRnarray() ]();
  }
#endif

  //--------------------------------------------------------------------------

  RandomNumberKernelBase::~RandomNumberKernelBase()
  {
#ifdef __CUDACC__
    // Deallocate the host buffer for the output random numbers
    if ( m_hstRnarray ) checkCuda( cudaFreeHost( m_hstRnarray ) );
    // Deallocate the device buffer for the output random numbers
    checkCuda( cudaFree( m_devRnarray ) );
#else
    // Deallocate the host buffer for the output random numbers
    ::operator delete( m_hstRnarray, std::align_val_t{ cppAlign } );
#endif
  }

  //--------------------------------------------------------------------------

  void RandomNumberKernelBase::seedGenerator( const int seed )
  {
    if ( seed <= 0 ) throw std::runtime_error( "RandomNumberKernelBase::seedGenerator expects a seed >0" );
    m_seed = seed;
  }

  //--------------------------------------------------------------------------

  void RandomNumberKernelBase::generateRnarray()
  {
    if ( m_seed <= 0 ) throw std::runtime_error( "RandomNumberKernelBase::generateRnarray expects a seed >0" );
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void RandomNumberKernelBase::copyHstRnarrayToDevRnarray()
  {
    const int nRnarray = np4 * nparf * m_nevt;
    const int nbytesRnarray = nRnarray * sizeof(fptype);
    if ( !m_hstRnarray )
      throw std::logic_error( "Cannot copyHsrRnarrayToDevRnarray unless random numbers are generated on the host" );
    checkCuda( cudaMemcpy( m_devRnarray, m_hstRnarray, nbytesRnarray, cudaMemcpyHostToDevice ) );
  }
#endif

  //--------------------------------------------------------------------------

  CommonRandomKernel::CommonRandomKernel( const int nevt )
    : RandomNumberKernelBase( nevt )
  {
  }

  //--------------------------------------------------------------------------

  CommonRandomKernel::~CommonRandomKernel()
  {
  }

  //--------------------------------------------------------------------------

  void CommonRandomKernel::generateRnarray()
  {
    RandomNumberKernelBase::generateRnarray(); // check if m_seed > 0
    std::vector<double> rnd = CommonRandomNumbers::generate<double>( nRnarray(), m_seed ); // NB: HARDCODED DOUBLE!
    std::copy( rnd.begin(), rnd.end(), m_hstRnarray ); // NB: this may imply a conversion from double to float
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HAS_NO_CURAND

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

  //--------------------------------------------------------------------------

  CurandRandomKernel::CurandRandomKernel( const int nevt, RandomNumberMode mode )
    : RandomNumberKernelBase( nevt )
    , m_mode( mode )
    , m_rnGen()
  {
    if ( m_mode != RandomNumberMode::CurandDevice && m_mode != RandomNumberMode::CurandHost )
      throw std::runtime_error( "CurandRandomKernel only supports CurandDevice and CurandHost" );
#ifndef __CUDACC__
    if ( m_mode == RandomNumberMode::CurandDevice )
      throw std::runtime_error( "CurandRandomKernel does not support CurandDevice on CPUs" );
#endif
    // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
    const curandRngType_t type = CURAND_RNG_PSEUDO_MTGP32;          // 0.00082s | 0.00064s (FOR FAST TESTS)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_XORWOW;        // 0.049s   | 0.0016s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MRG32K3A;      // 0.71s    | 0.0012s  (better but slower, especially in c++)
    //const curandRngType_t type = CURAND_RNG_PSEUDO_MT19937;       // 21s      | 0.021s
    //const curandRngType_t type = CURAND_RNG_PSEUDO_PHILOX4_32_10; // 0.024s   | 0.00026s (used to segfault?)
    if ( m_mode == RandomNumberMode::CurandDevice )
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

  CurandRandomKernel::~CurandRandomKernel()
  {
    checkCurand( curandDestroyGenerator( m_rnGen ) );
  }

  //--------------------------------------------------------------------------

  void CurandRandomKernel::seedGenerator( const int seed )
  {
    RandomNumberKernelBase::seedGenerator( seed ); // check if seed > 0 and set m_seed
    //printf( "seedGenerator: seed %lld\n", seed );
    checkCurand( curandSetPseudoRandomGeneratorSeed( m_rnGen, seed ) );
  }

  //--------------------------------------------------------------------------

  void CurandRandomKernel::generateRnarray()
  {
    RandomNumberKernelBase::generateRnarray(); // check if m_seed > 0
    fptype* outRnarray = nullptr;
#ifdef __CUDACC__
    if ( m_mode == RandomNumberMode::CurandDevice )
      outRnarray = m_devRnarray;
    else outRnarray = m_hstRnarray;
#else
    if ( m_mode == RandomNumberMode::CurandDevice )
      throw std::runtime_error( "CurandRandomKernel does not support CurandDevice on CPUs" );
    else outRnarray = m_hstRnarray;
#endif
#if defined MGONGPU_FPTYPE_DOUBLE
    checkCurand( curandGenerateUniformDouble( m_rnGen, outRnarray, nRnarray() ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkCurand( curandGenerateUniform( m_rnGen, outRnarray, nRnarray() ) );
#endif
    //for ( int i=0; i<8; i++ ) printf("%f %f %f %f\n",outRrnarray[i*4],outRrnarray[i*4+2],outRrnarray[i*4+2],outRrnarray[i*4+3]);
  }

  //--------------------------------------------------------------------------
#endif

}
