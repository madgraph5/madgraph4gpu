// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "MemoryBuffers.h"
#include "RandomNumberKernels.h"

#include <cassert>

#ifndef MGONGPU_HAS_NO_ROCRAND /* clang-format off */
#include "rocrand.h"
#define checkRocRand( code ){ assertRocRand( code, __FILE__, __LINE__ ); }
inline void assertRocRand( rocrandStatus_t code, const char *file, int line, bool abort = true )
{
  if ( code != ROCRAND_STATUS_SUCCESS )
  {
    printf( "RocRandAssert: %s:%d code=%d\n", file, line, code );
    if ( abort ) assert( code == ROCRAND_STATUS_SUCCESS );
  }
}
#endif /* clang-format on */

#ifdef MGONGPUCPP_HIPCC
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------
#ifndef MGONGPU_HAS_NO_ROCRAND
  RocRandRandomNumberKernel::RocRandRandomNumberKernel( BufferRndNumMomenta& rnarray, const bool onDevice )
    : RandomNumberKernelBase( rnarray )
    , m_isOnDevice( onDevice )
  {
    if( m_isOnDevice )
    {
#ifdef MGONGPUCPP_HIPCC
      if( !m_rnarray.isOnDevice() )
        throw std::runtime_error( "RocRandRandomNumberKernel on device with a host random number array" );
#else
      throw std::runtime_error( "RocRandRandomNumberKernel does not support RocRandDevice on CPU host" );
#endif
    }
    else
    {
      if( m_rnarray.isOnDevice() )
        throw std::runtime_error( "RocRandRandomNumberKernel on host with a device random number array" );
    }
    createGenerator();
  }

  //--------------------------------------------------------------------------

  RocRandRandomNumberKernel::~RocRandRandomNumberKernel()
  {
    destroyGenerator();
  }

  //--------------------------------------------------------------------------

  void RocRandRandomNumberKernel::seedGenerator( const unsigned int seed )
  {
    if( m_isOnDevice )
    {
      destroyGenerator(); // workaround for #429
      createGenerator();  // workaround for #429
    }
    //printf( "seedGenerator: seed %d\n", seed );
    checkRocRand( rocrand_set_seed( m_rnGen, seed ) );
  }

  //--------------------------------------------------------------------------

  void RocRandRandomNumberKernel::createGenerator()
  {
    // [NB Timings are for GenRnGen host|device (cpp|cuda) generation of 256*32*1 events with nproc=1: rn(0) is host=0.0012s]
    const rocrand_rng_type type = ROCRAND_RNG_PSEUDO_MTGP32; // 0.00082s | 0.00064s (FOR FAST TESTS) NO TESTS DONE WITH ROCRAND, BUT THIS IS FASTEST IN CUDA
    //const rocrand_rng_type type = ROCRAND_RNG_PSEUDO_XORWOW;
    //const rocrand_rng_type type = ROCRAND_RNG_PSEUDO_MRG32K3A;
    //const rocrand_rng_type type = ROCRAND_RNG_PSEUDO_MT19937;
    //const rocrand_rng_type type = ROCRAND_RNG_PSEUDO_PHILOX4_32_10;
    if( m_isOnDevice )
    {
      checkRocRand( rocrand_create_generator( &m_rnGen, type ) );
    }
    else
    {
      checkRocRand( rocrand_create_generator_host( &m_rnGen, type ) );
    }
    // No RocRAND equivalent for curandSetGeneratorOrdering
  }

  //--------------------------------------------------------------------------

  void RocRandRandomNumberKernel::destroyGenerator()
  {
    checkRocRand( rocrand_destroy_generator( m_rnGen ) );
  }

  //--------------------------------------------------------------------------

  void RocRandRandomNumberKernel::generateRnarray()
  {
#if defined MGONGPU_FPTYPE_DOUBLE
    checkRocRand( rocrand_generate_uniform_double( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#elif defined MGONGPU_FPTYPE_FLOAT
    checkRocRand( rocrand_generate_uniform( m_rnGen, m_rnarray.data(), m_rnarray.size() ) );
#endif
    /* 
    printf( "\nRocRandRandomNumberKernel::generateRnarray size = %d\n", (int)m_rnarray.size() );
    fptype* data = m_rnarray.data();
#ifdef MGONGPUCPP_GPUIMPL
    if( m_rnarray.isOnDevice() )
    {
      data = new fptype[m_rnarray.size()]();
      checkRoc( rocMemcpy( data, m_rnarray.data(), m_rnarray.bytes(), rocMemcpyDeviceToHost ) );
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