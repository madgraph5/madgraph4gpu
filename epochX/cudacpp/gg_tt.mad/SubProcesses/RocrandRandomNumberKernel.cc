// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "MemoryBuffers.h"
#include "RandomNumberKernels.h"

#include <cassert>

#ifndef MGONGPU_HAS_NO_ROCRAND /* clang-format off */
#include <rocrand/rocrand.hpp>
#define checkRocrand( code ){ assertRocrand( code, __FILE__, __LINE__ ); }
inline void assertRocrand( rocrandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != ROCRAND_STATUS_SUCCESS )
  {
    printf( "RocrandAssert: %s:%d code=%d\n", file, line, code );
    if ( abort ) assert( code == ROCRAND_STATUS_SUCCESS );
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
      throw std::runtime_error( "RocrandRandomNumberKernel does not support RocrandDevice on CPU host" );
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
    checkRocrand( rocrandSetPseudoRandomGeneratorSeed( m_rnGen, seed ) );
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::createGenerator()
  {
    // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
    const rocrandRngType_t type = ROCRAND_RNG_PSEUDO_MTGP32; //          0.00082s | 0.00064s (FOR FAST TESTS)
    //const rocrandRngType_t type = ROCRAND_RNG_PSEUDO_XORWOW;        // 0.049s   | 0.0016s
    //const rocrandRngType_t type = ROCRAND_RNG_PSEUDO_MRG32K3A;      // 0.71s    | 0.0012s  (better but slower, especially in c++)
    //const rocrandRngType_t type = ROCRAND_RNG_PSEUDO_MT19937;       // 21s      | 0.021s
    //const rocrandRngType_t type = ROCRAND_RNG_PSEUDO_PHILOX4_32_10; // 0.024s   | 0.00026s (used to segfault?)
    if( m_isOnDevice )
    {
      checkRocrand( rocrandCreateGenerator( &m_rnGen, type ) );
    }
    else
    {
      checkRocrand( rocrandCreateGeneratorHost( &m_rnGen, type ) );
    }
    //checkRocrand( rocrandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_LEGACY ) ); // fails with code=104 (see #429)
    checkRocrand( rocrandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_BEST ) );
    //checkRocrand( rocrandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_DYNAMIC ) ); // fails with code=104 (see #429)
    //checkRocrand( rocrandSetGeneratorOrdering( *&m_rnGen, ROCRAND_ORDERING_PSEUDO_SEEDED ) ); // fails with code=104 (see #429)
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::destroyGenerator()
  {
    checkRocrand( rocrandDestroyGenerator( m_rnGen ) );
  }

  //--------------------------------------------------------------------------

  void RocrandRandomNumberKernel::generateRnarray()
  {
#if defined MGONGPU_FPTYPE_DOUBLE
    checkRocrand( rocrandGenerateUniformDouble( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkRocrand( rocrandGenerateUniform( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
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
