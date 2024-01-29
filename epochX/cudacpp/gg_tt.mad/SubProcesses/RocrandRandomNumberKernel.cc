// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "MemoryBuffers.h"
#include "RandomNumberKernels.h"

#include <cassert>

#ifndef MGONGPU_HAS_NO_ROCRAND /* clang-format off */
//#include <hiprand/hiprand.hpp>
#include "hiprand.h"
#define checkHiprand( code ){ assertHiprand( code, __FILE__, __LINE__ ); }
inline void assertHiprand( hiprandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != HIPRAND_STATUS_SUCCESS )
  {
    printf( "HiprandAssert: %s:%d code=%d\n", file, line, code );
    if ( abort ) assert( code == HIPRAND_STATUS_SUCCESS );
  }
}
#endif /* clang-format on */

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------
#ifndef MGONGPU_HAS_NO_ROCRAND
  RocrandRandomNumberKernel::RocrandRandomNumberKernel( BufferRndNumMomenta& rnarray, const bool onDevice )
    : RandomNumberKernelBase( rnarray )
    , m_isOnDevice( onDevice )
  {
    if( m_isOnDevice )
    {
#ifdef MGONGPUCPP_GPUIMPL
      if( !m_rnarray.isOnDevice() )
        throw std::runtime_error( "RocrandRandomNumberKernel on device with a host random number array" );
#else
      throw std::runtime_error( "RocrandRandomNumberKernel does not support HiprandDevice on CPU host" );
#endif
    }
    else
    {
      if( m_rnarray.isOnDevice() )
        throw std::runtime_error( "RocrandRandomNumberKernel on host with a device random number array" );
    }
    createGenerator();
  }

  //--------------------------------------------------------------------------

  RocrandRandomNumberKernel::~RocrandRandomNumberKernel()
  {
    destroyGenerator();
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::seedGenerator( const unsigned int seed )
  {
    if( m_isOnDevice )
    {
      destroyGenerator(); // workaround for #429
      createGenerator();  // workaround for #429
    }
    //printf( "seedGenerator: seed %d\n", seed );
    checkHiprand( hiprandSetPseudoRandomGeneratorSeed( m_rnGen, seed ) );
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::createGenerator()
  {
    // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
    const hiprandRngType_t type = HIPRAND_RNG_PSEUDO_MTGP32; //          0.00082s | 0.00064s (FOR FAST TESTS)
    //const hiprandRngType_t type = HIPRAND_RNG_PSEUDO_XORWOW;        // 0.049s   | 0.0016s
    //const hiprandRngType_t type = HIPRAND_RNG_PSEUDO_MRG32K3A;      // 0.71s    | 0.0012s  (better but slower, especially in c++)
    //const hiprandRngType_t type = HIPRAND_RNG_PSEUDO_MT19937;       // 21s      | 0.021s
    //const hiprandRngType_t type = HIPRAND_RNG_PSEUDO_PHILOX4_32_10; // 0.024s   | 0.00026s (used to segfault?)
    if( m_isOnDevice )
    {
      checkHiprand( hiprandCreateGenerator( &m_rnGen, type ) );
    }
    else
    {
      checkHiprand( hiprandCreateGeneratorHost( &m_rnGen, type ) );
    }
    /*
    // FIXME: implement hiprand/rocrand ordering...
    //checkHiprand( hiprandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_LEGACY ) ); // fails in curand (see #429)
    checkHiprand( hiprandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_BEST ) );
    //checkHiprand( hiprandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_DYNAMIC ) ); // fails in curand (see #429)
    //checkHiprand( hiprandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_SEEDED ) ); // fails in curand (see #429)
    */
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::destroyGenerator()
  {
    checkHiprand( hiprandDestroyGenerator( m_rnGen ) );
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::generateRnarray()
  {
#if defined MGONGPU_FPTYPE_DOUBLE
    checkHiprand( hiprandGenerateUniformDouble( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkHiprand( hiprandGenerateUniform( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#endif
    /*
    printf( "\nRocrandRandomNumberKernel::generateRnarray size = %d\n", (int)m_rnarray.size() );
    fptype* data = m_rnarray.data();
#ifdef MGONGPUCPP_GPUIMPL
    if( m_rnarray.isOnDevice() )
    {
      data = new fptype[m_rnarray.size()]();
      checkCuda( cudaMemcpy( data, m_rnarray.data(), m_rnarray.bytes(), cudaMemcpyDeviceToHost ) );
    }
#endif
    for( int i = 0; i < ( (int)m_rnarray.size() / 4 ); i++ )
      printf( "[%4d] %f %f %f %f\n", i * 4, data[i * 4], data[i * 4 + 2], data[i * 4 + 2], data[i * 4 + 3] );
#ifdef MGONGPUCPP_GPUIMPL
    if( m_rnarray.isOnDevice() ) delete[] data;
#endif
    */
  }

  //--------------------------------------------------------------------------
#endif
}
