// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Hageboeck (Nov 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
//----------------------------------------------------------------------------
// Use ./runTest.exe --gtest_filter=*xxx to run only testxxx.cc tests
//----------------------------------------------------------------------------

#include "mgOnGpuConfig.h"

#include "CPPProcess.h"
#include "MadgraphTest.h"
#include "MatrixElementKernels.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryBuffers.h"
#include "RamboSamplingKernels.h"
#include "RandomNumberKernels.h"
#include "epoch_process_id.h"

#ifdef __CUDACC__
using namespace mg5amcGpu;
#else
using namespace mg5amcCpu;
#endif

struct CUDA_CPU_TestBase : public TestDriverBase
{
  static constexpr int neppM = MemoryAccessMomenta::neppM; // AOSOA layout
  static constexpr int np4 = CPPProcess::np4;
  static constexpr int npar = CPPProcess::npar;
  static_assert( gputhreads % neppM == 0, "ERROR! #threads/block should be a multiple of neppM" );
  static_assert( gputhreads <= mgOnGpu::ntpbMAX, "ERROR! #threads/block should be <= ntpbMAX" );
  CUDA_CPU_TestBase( const std::string& refFileName )
    : TestDriverBase( npar, refFileName ) {}
};

#ifndef __CUDACC__
struct CPUTest : public CUDA_CPU_TestBase
{
  // Struct data members (process, and memory structures for random numbers, momenta, matrix elements and weights on host and device)
  // [NB the hst/dev memory arrays must be initialised in the constructor, see issue #290]
  CPPProcess process;
  HostBufferRndNumMomenta hstRndmom;
  HostBufferMomenta hstMomenta;
  HostBufferGs hstGs;
  HostBufferRndNumHelicity hstRndHel;
  HostBufferRndNumColor hstRndCol;
  HostBufferWeights hstWeights;
  HostBufferMatrixElements hstMatrixElements;
  HostBufferSelectedHelicity hstSelHel;
  HostBufferSelectedColor hstSelCol;
  HostBufferHelicityMask hstIsGoodHel;

  // Create a process object
  // Read param_card and set parameters
  // ** WARNING EVIL EVIL **
  // The CPPProcess constructor has side effects on the globals Proc::cHel, which is needed in ME calculations.
  // Don't remove!
  CPUTest( const std::string& refFileName )
    : CUDA_CPU_TestBase( refFileName )
    , process( /*verbose=*/false )
    , hstRndmom( nevt )
    , hstMomenta( nevt )
    , hstGs( nevt )
    , hstRndHel( nevt )
    , hstRndCol( nevt )
    , hstWeights( nevt )
    , hstMatrixElements( nevt )
    , hstSelHel( nevt )
    , hstSelCol( nevt )
    , hstIsGoodHel( CPPProcess::ncomb )
  {
    // FIXME: the process instance can happily go out of scope because it is only needed to read parameters?
    // FIXME: the CPPProcess should really be a singleton?
    process.initProc( "../../Cards/param_card.dat" );
  }

  virtual ~CPUTest() {}

  void prepareRandomNumbers( unsigned int iiter ) override
  {
    CommonRandomNumberKernel rnk( hstRndmom );
    rnk.seedGenerator( 1337 + iiter );
    rnk.generateRnarray();
  }

  void prepareMomenta( fptype energy ) override
  {
    RamboSamplingKernelHost rsk( energy, hstRndmom, hstMomenta, hstWeights, nevt );
    // --- 2a. Fill in momenta of initial state particles on the device
    rsk.getMomentaInitial();
    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    rsk.getMomentaFinal();
  }

  void runSigmaKin( std::size_t iiter ) override
  {
    constexpr fptype fixedG = 1.2177157847767195; // fixed G for aS=0.118 (hardcoded for now in check_sa.cc, fcheck_sa.f, runTest.cc)
    for( unsigned int i = 0; i < nevt; ++i ) hstGs[i] = fixedG;
    MatrixElementKernelHost mek( hstMomenta, hstGs, hstRndHel, hstRndCol, hstMatrixElements, hstSelHel, hstSelCol, nevt );
    if( iiter == 0 ) mek.computeGoodHelicities();
    constexpr unsigned int channelId = 0; // TEMPORARY? disable multi-channel in runTest.exe #466
    mek.computeMatrixElements( channelId );
  }

  fptype getMomentum( std::size_t ievt, unsigned int ipar, unsigned int ip4 ) const override
  {
    assert( ipar < npar );
    assert( ip4 < np4 );
    return MemoryAccessMomenta::ieventAccessIp4IparConst( hstMomenta.data(), ievt, ip4, ipar );
  }

  fptype getMatrixElement( std::size_t ievt ) const override
  {
    return MemoryAccessMatrixElements::ieventAccessConst( hstMatrixElements.data(), ievt );
  }
};
#endif

#ifdef __CUDACC__
struct CUDATest : public CUDA_CPU_TestBase
{
  // Reset the device when our test goes out of scope. Note that this should happen after
  // the frees, i.e. be declared before the pointers to device memory.
  struct DeviceReset
  {
    ~DeviceReset()
    {
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
  } deviceResetter;

  // Struct data members (process, and memory structures for random numbers, momenta, matrix elements and weights on host and device)
  // [NB the hst/dev memory arrays must be initialised in the constructor, see issue #290]
  CPPProcess process;
  PinnedHostBufferRndNumMomenta hstRndmom;
  PinnedHostBufferMomenta hstMomenta;
  PinnedHostBufferGs hstGs;
  PinnedHostBufferRndNumHelicity hstRndHel;
  PinnedHostBufferRndNumColor hstRndCol;
  PinnedHostBufferWeights hstWeights;
  PinnedHostBufferMatrixElements hstMatrixElements;
  PinnedHostBufferSelectedHelicity hstSelHel;
  PinnedHostBufferSelectedColor hstSelCol;
  PinnedHostBufferHelicityMask hstIsGoodHel;
  DeviceBufferRndNumMomenta devRndmom;
  DeviceBufferMomenta devMomenta;
  DeviceBufferGs devGs;
  DeviceBufferRndNumHelicity devRndHel;
  DeviceBufferRndNumColor devRndCol;
  DeviceBufferWeights devWeights;
  DeviceBufferMatrixElements devMatrixElements;
  DeviceBufferSelectedHelicity devSelHel;
  DeviceBufferSelectedColor devSelCol;
  DeviceBufferHelicityMask devIsGoodHel;

  // Create a process object
  // Read param_card and set parameters
  // ** WARNING EVIL EVIL **
  // The CPPProcess constructor has side effects on the globals Proc::cHel, which is needed in ME calculations.
  // Don't remove!
  CUDATest( const std::string& refFileName )
    : CUDA_CPU_TestBase( refFileName )
    , process( /*verbose=*/false )
    , hstRndmom( nevt )
    , hstMomenta( nevt )
    , hstGs( nevt )
    , hstRndHel( nevt )
    , hstRndCol( nevt )
    , hstWeights( nevt )
    , hstMatrixElements( nevt )
    , hstSelHel( nevt )
    , hstSelCol( nevt )
    , hstIsGoodHel( CPPProcess::ncomb )
    , devRndmom( nevt )
    , devMomenta( nevt )
    , devGs( nevt )
    , devRndHel( nevt )
    , devRndCol( nevt )
    , devWeights( nevt )
    , devMatrixElements( nevt )
    , devSelHel( nevt )
    , devSelCol( nevt )
    , devIsGoodHel( CPPProcess::ncomb )
  {
    // FIXME: the process instance can happily go out of scope because it is only needed to read parameters?
    // FIXME: the CPPProcess should really be a singleton?
    process.initProc( "../../Cards/param_card.dat" );
  }

  virtual ~CUDATest() {}

  void prepareRandomNumbers( unsigned int iiter ) override
  {
    CommonRandomNumberKernel rnk( hstRndmom );
    rnk.seedGenerator( 1337 + iiter );
    rnk.generateRnarray();
    copyDeviceFromHost( devRndmom, hstRndmom );
  }

  void prepareMomenta( fptype energy ) override
  {
    RamboSamplingKernelDevice rsk( energy, devRndmom, devMomenta, devWeights, gpublocks, gputhreads );
    // --- 2a. Fill in momenta of initial state particles on the device
    rsk.getMomentaInitial();
    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    rsk.getMomentaFinal();
    // --- 2c. CopyDToH Weights
    copyHostFromDevice( hstWeights, devWeights );
    // --- 2d. CopyDToH Momenta
    copyHostFromDevice( hstMomenta, devMomenta );
  }

  void runSigmaKin( std::size_t iiter ) override
  {
    constexpr fptype fixedG = 1.2177157847767195; // fixed G for aS=0.118 (hardcoded for now in check_sa.cc, fcheck_sa.f, runTest.cc)
    for( unsigned int i = 0; i < nevt; ++i ) hstGs[i] = fixedG;
    copyDeviceFromHost( devGs, hstGs ); // BUG FIX #566
    MatrixElementKernelDevice mek( devMomenta, devGs, devRndHel, devRndCol, devMatrixElements, devSelHel, devSelCol, gpublocks, gputhreads );
    if( iiter == 0 ) mek.computeGoodHelicities();
    constexpr unsigned int channelId = 0; // TEMPORARY? disable multi-channel in runTest.exe #466
    mek.computeMatrixElements( channelId );
    copyHostFromDevice( hstMatrixElements, devMatrixElements );
  }

  fptype getMomentum( std::size_t ievt, unsigned int ipar, unsigned int ip4 ) const override
  {
    assert( ipar < npar );
    assert( ip4 < np4 );
    return MemoryAccessMomenta::ieventAccessIp4IparConst( hstMomenta.data(), ievt, ip4, ipar );
  }

  fptype getMatrixElement( std::size_t ievt ) const override
  {
    return MemoryAccessMatrixElements::ieventAccessConst( hstMatrixElements.data(), ievt );
  }
};
#endif

// Use two levels of macros to force stringification at the right level
// (see https://gcc.gnu.org/onlinedocs/gcc-3.0.1/cpp_3.html#SEC17 and https://stackoverflow.com/a/3419392)
// Google macro is in https://github.com/google/googletest/blob/master/googletest/include/gtest/gtest-param-test.h
#define TESTID_CPU( s ) s##_CPU
#define XTESTID_CPU( s ) TESTID_CPU( s )
#define MG_INSTANTIATE_TEST_SUITE_CPU( prefix, test_suite_name ) \
INSTANTIATE_TEST_SUITE_P( prefix, \
                          test_suite_name, \
                          testing::Values( new CPUTest( MG_EPOCH_REFERENCE_FILE_NAME ) ) );
#define TESTID_GPU( s ) s##_GPU
#define XTESTID_GPU( s ) TESTID_GPU( s )
#define MG_INSTANTIATE_TEST_SUITE_GPU( prefix, test_suite_name ) \
INSTANTIATE_TEST_SUITE_P( prefix, \
                          test_suite_name, \
                          testing::Values( new CUDATest( MG_EPOCH_REFERENCE_FILE_NAME ) ) );

#ifdef __CUDACC__
MG_INSTANTIATE_TEST_SUITE_GPU( XTESTID_GPU( MG_EPOCH_PROCESS_ID ), MadgraphTest );
#else
MG_INSTANTIATE_TEST_SUITE_CPU( XTESTID_CPU( MG_EPOCH_PROCESS_ID ), MadgraphTest );
#endif
