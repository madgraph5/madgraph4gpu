#include "RandomNumberKernel.h"

#include "checkCuda.h"
#include "CommonRandomNumbers.h"
#include "rambo.h"

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

#ifndef MGONGPU_COMMONRAND_ONHOST
  CurandRandomKernel::CurandRandomKernel( const int nevt, RandomNumberMode mode )
    : RandomNumberKernelBase( nevt )
    , m_mode( mode )
    , m_rnGen()
  {
#ifdef __CUDACC__
    grambo2toNm0::createGenerator( &m_rnGen );
#else
    if ( m_mode == RandomNumberMode::CurandDevice )
      throw std::runtime_error( "CurandRandomKernel does not support CurandDevice on CPUs" );
    rambo2toNm0::createGenerator( &m_rnGen );
#endif
  }

  //--------------------------------------------------------------------------

  CurandRandomKernel::~CurandRandomKernel()
  {
#ifdef __CUDACC__
    grambo2toNm0::destroyGenerator( m_rnGen );
#else
    rambo2toNm0::destroyGenerator( m_rnGen );
#endif
  }

  //--------------------------------------------------------------------------

  void CurandRandomKernel::seedGenerator( const int seed )
  {
    RandomNumberKernelBase::seedGenerator( seed ); // check if seed > 0 and set m_seed
#ifdef __CUDACC__
    grambo2toNm0::seedGenerator( m_rnGen, seed );
#else
    rambo2toNm0::seedGenerator( m_rnGen, seed );
#endif
  }

  //--------------------------------------------------------------------------

  void CurandRandomKernel::generateRnarray()
  {
    RandomNumberKernelBase::generateRnarray(); // check if m_seed > 0
#ifdef __CUDACC__
    if ( m_mode == RandomNumberMode::CurandDevice )
      grambo2toNm0::generateRnarray( m_rnGen, m_devRnarray, m_nevt );
    else
      grambo2toNm0::generateRnarray( m_rnGen, m_hstRnarray, m_nevt );
#else
    if ( m_mode == RandomNumberMode::CurandDevice )
      throw std::runtime_error( "CurandRandomKernel does not support CurandDevice on CPUs" );
    else
      rambo2toNm0::generateRnarray( m_rnGen, m_hstRnarray, m_nevt );
#endif
  }

  //--------------------------------------------------------------------------
#endif

}
