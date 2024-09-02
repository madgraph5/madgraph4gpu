// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Hageboeck (Nov 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, J. Teig, A. Valassi (2020-2024) for the MG5aMC CUDACPP plugin.
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

#ifdef MGONGPUCPP_GPUIMPL
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

#ifndef MGONGPUCPP_GPUIMPL
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
  std::unique_ptr<MatrixElementKernelBase> pmek;

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
    , pmek( new MatrixElementKernelHost( hstMomenta, hstGs, hstRndHel, hstRndCol, hstMatrixElements, hstSelHel, hstSelCol, nevt ) )
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
    if( iiter == 0 ) pmek->computeGoodHelicities();
    constexpr unsigned int channelId = 0; // TEMPORARY? disable multi-channel in runTest.exe #466
    pmek->computeMatrixElements( channelId );
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

#ifdef MGONGPUCPP_GPUIMPL
struct CUDATest : public CUDA_CPU_TestBase
{
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
  std::unique_ptr<MatrixElementKernelBase> pmek;

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
    , pmek( new MatrixElementKernelDevice( devMomenta, devGs, devRndHel, devRndCol, devMatrixElements, devSelHel, devSelCol, gpublocks, gputhreads ) )
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
    if( iiter == 0 ) pmek->computeGoodHelicities();
    constexpr unsigned int channelId = 0; // TEMPORARY? disable multi-channel in runTest.exe #466
    pmek->computeMatrixElements( channelId );
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
#endif /* clang-format off */

// AV July 2024 much simpler class structure without the presently-unnecessary googletest templates
// This is meant as a workaround to prevent not-understood segfault #907 when adding a second test
#ifdef MGONGPUCPP_GPUIMPL
// CUDA test 1
CUDATest cudaDriver1( MG_EPOCH_REFERENCE_FILE_NAME );
MadgraphTest mgTest1( cudaDriver1 );
#define TESTID1( s ) s##_GPU_MADGRAPH1
#define XTESTID1( s ) TESTID1( s )
// CUDA test 2
//CUDATest cudaDriver2( MG_EPOCH_REFERENCE_FILE_NAME );
//MadgraphTest mgTest2( cudaDriver2 );
//#define TESTID2( s ) s##_GPU_MADGRAPH2
//#define XTESTID2( s ) TESTID2( s )
#else
// CPU test 1
CPUTest cppDriver1( MG_EPOCH_REFERENCE_FILE_NAME );
MadgraphTest mgTest1( cppDriver1 );
#define TESTID1( s ) s##_CPU_MADGRAPH1
#define XTESTID1( s ) TESTID1( s )
// CPU test 2
//CPUTest cppDriver2( MG_EPOCH_REFERENCE_FILE_NAME );
//MadgraphTest mgTest2( cppDriver2 );
//#define TESTID2( s ) s##_CPU_MADGRAPH2
//#define XTESTID2( s ) TESTID2( s )
#endif
// Instantiate Google test 1
TEST( XTESTID1( MG_EPOCH_PROCESS_ID ), compareMomAndME )
{
  mgTest1.CompareMomentaAndME( *this );
}
// Instantiate Google test 2
//TEST( XTESTID2( MG_EPOCH_PROCESS_ID ), compareMomAndME )
//{
//  mgTest2.CompareMomentaAndME( *this );
//}
/* clang-format on */
