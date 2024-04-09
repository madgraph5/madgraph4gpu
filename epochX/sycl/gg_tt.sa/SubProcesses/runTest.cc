// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Hageboeck (Nov 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, O. Mattelaer, S. Roiser, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
//
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.

#include "epoch_process_id.h"

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "MadgraphTest.h"

#include "CommonRandomNumbers.h"
#include "CPPProcess.h"
#include "Memory.h"
#ifdef __CUDACC__
#include "rambo.cc"
#else
#include "rambo.h"
#endif

#ifdef __CUDACC__
template<typename T = fptype>
using unique_ptr_dev = std::unique_ptr<T, CudaDevDeleter<T>>;
template<typename T = fptype>
using unique_ptr_host = std::unique_ptr<T[], CudaHstDeleter<T>>;
#else
template<typename T = fptype>
using unique_ptr_host = std::unique_ptr<T[]>;
#endif

struct CUDA_CPU_TestBase : public TestDriverBase<fptype> {

  static_assert( gputhreads%mgOnGpu::neppR == 0, "ERROR! #threads/block should be a multiple of neppR" );
  static_assert( gputhreads%mgOnGpu::neppM == 0, "ERROR! #threads/block should be a multiple of neppM" );
  static_assert( gputhreads <= mgOnGpu::ntpbMAX, "ERROR! #threads/block should be <= ntpbMAX" );

  const std::size_t nRnarray{ CPPPROCESS_NP4*CPPPROCESS_NPARF*nevt }; // AOSOA layout with nevt=npagR*neppR events per iteration
  const std::size_t nMomenta{ CPPPROCESS_NP4*CPPPROCESS_NPAR*nevt }; // AOSOA layout with nevt=npagM*neppM events per iteration
  const std::size_t nWeights{ nevt };
  const std::size_t nMEs    { nevt };

  CUDA_CPU_TestBase() :
    TestDriverBase()
  {
    TestDriverBase::nparticle = CPPPROCESS_NPAR;
  }

};

#ifndef __CUDACC__
struct CPUTest : public CUDA_CPU_TestBase {

  Proc::CPPProcess process;

  // --- 0b. Allocate memory structures
  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  unique_ptr_host<fptype   > hstRnarray  { hstMakeUnique<fptype   >( nRnarray ) }; // AOSOA[npagR][CPPPROCESS_NPARF][CPPPROCESS_NP4][neppR]
  unique_ptr_host<fptype_sv> hstMomenta  { hstMakeUnique<fptype_sv>( nMomenta ) }; // AOSOA[npagM][CPPPROCESS_NPAR][CPPPROCESS_NP4][neppM]
  unique_ptr_host<bool     > hstIsGoodHel{ hstMakeUnique<bool     >( CPPPROCESS_NCOMB ) };
  unique_ptr_host<fptype   > hstWeights  { hstMakeUnique<fptype   >( nWeights ) };
  unique_ptr_host<fptype_sv> hstMEs      { hstMakeUnique<fptype_sv>( nMEs ) }; // AOSOA[npagM][neppM]

  // Create a process object
  // Read param_card and set parameters
  // ** WARNING EVIL EVIL **
  // The CPPProcess constructor has side effects on the globals Proc::cHel, which is needed in ME calculations.
  // Don't remove!
  CPUTest() :
    CUDA_CPU_TestBase(),
    process(niter, gpublocks, gputhreads, /*verbose=*/false)
  {
    process.initProc("../../Cards/param_card.dat");
  }
  virtual ~CPUTest() { }

  void prepareRandomNumbers(unsigned int iiter) override {
    std::vector<double> rnd = CommonRandomNumbers::generate<double>(nRnarray, 1337 + iiter); // NB: HARDCODED DOUBLE!
    std::copy(rnd.begin(), rnd.end(), hstRnarray.get()); // NB: this may imply a conversion from double to float
  }

  void prepareMomenta(fptype energy) override {
    // --- 2a. Fill in momenta of initial state particles on the device
    rambo2toNm0::getMomentaInitial( energy, hstMomenta.get(), nevt );
    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    rambo2toNm0::getMomentaFinal( energy, hstRnarray.get(), hstMomenta.get(), hstWeights.get(), nevt );
  }

  void runSigmaKin(std::size_t iiter) override {
    // --- 0d. SGoodHel
    if ( iiter == 0 )
    {
      // ... 0d1. Compute good helicity mask on the host
      Proc::sigmaKin_getGoodHel(hstMomenta.get(), hstMEs.get(), hstIsGoodHel.get(), nevt);
      // ... 0d2. Copy back good helicity list to static memory on the host
      Proc::sigmaKin_setGoodHel(hstIsGoodHel.get());
    }

    // --- 3a. SigmaKin
    Proc::sigmaKin(hstMomenta.get(), hstMEs.get(), nevt);
  }

  fptype getMomentum(std::size_t evtNo, unsigned int particle, unsigned int component) const override {
    using mgOnGpu::neppM;
    assert(component < CPPPROCESS_NP4);
    assert(particle  < CPPPROCESS_NPAR);
    const auto ipagM = evtNo / neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % neppM; // #event in the current eventpage in this iteration
#ifndef MGONGPU_CPPSIMD
    return hstMomenta[ipagM*CPPPROCESS_NPAR*CPPPROCESS_NP4*neppM + particle*CPPPROCESS_NP4*neppM + component*neppM + ieppM];
#else
    return hstMomenta[ipagM*CPPPROCESS_NPAR*CPPPROCESS_NP4 + particle*CPPPROCESS_NP4 + component][ieppM];
#endif
  };

  fptype getMatrixElement(std::size_t ievt) const override {
#ifndef MGONGPU_CPPSIMD
    return hstMEs[ievt];
#else
    return hstMEs[ievt/neppV][ievt%neppV];
#endif
  }

};
#endif

#ifdef __CUDACC__
struct CUDATest : public CUDA_CPU_TestBase {

  // Reset the device when our test goes out of scope. Note that this should happen after
  // the frees, i.e. be declared before the pointers to device memory.
  struct DeviceReset {
    ~DeviceReset() {
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
  } deviceResetter;

  // --- 0b. Allocate memory structures
  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  unique_ptr_host<fptype> hstRnarray  { hstMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][CPPPROCESS_NPARF][CPPPROCESS_NP4][neppR] (nevt=npagR*neppR)
  unique_ptr_host<fptype> hstMomenta  { hstMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][CPPPROCESS_NPAR][CPPPROCESS_NP4][neppM] (nevt=npagM*neppM)
  unique_ptr_host<bool  > hstIsGoodHel{ hstMakeUnique<bool  >( CPPPROCESS_NCOMB ) };
  unique_ptr_host<fptype> hstWeights  { hstMakeUnique<fptype>( nWeights ) };
  unique_ptr_host<fptype> hstMEs      { hstMakeUnique<fptype>( nMEs ) };

  unique_ptr_dev<fptype> devRnarray  { devMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][CPPPROCESS_NPARF][CPPPROCESS_NP4][neppR] (nevt=npagR*neppR)
  unique_ptr_dev<fptype> devMomenta  { devMakeUnique<fptype>( nMomenta ) };
  unique_ptr_dev<bool  > devIsGoodHel{ devMakeUnique<bool  >( CPPPROCESS_NCOMB ) };
  unique_ptr_dev<fptype> devWeights  { devMakeUnique<fptype>( nWeights ) };
  unique_ptr_dev<fptype> devMEs      { devMakeUnique<fptype>( nMEs )     };

  gProc::CPPProcess process;

  // Create a process object
  // Read param_card and set parameters
  // ** WARNING EVIL EVIL **
  // The CPPProcess constructor has side effects on the globals Proc::cHel, which is needed in ME calculations.
  // Don't remove!
  CUDATest() :
    CUDA_CPU_TestBase(),
    process(niter, gpublocks, gputhreads, /*verbose=*/false)
  {
    process.initProc("../../Cards/param_card.dat");
  }

  virtual ~CUDATest() { }

  void prepareRandomNumbers(unsigned int iiter) override {
    std::vector<double> rnd = CommonRandomNumbers::generate<double>(nRnarray, 1337 + iiter); // NB: HARDCODED DOUBLE!
    std::copy(rnd.begin(), rnd.end(), hstRnarray.get()); // NB: this may imply a conversion from double to float
    checkCuda( cudaMemcpy( devRnarray.get(), hstRnarray.get(),
                           nRnarray * sizeof(decltype(devRnarray)::element_type), cudaMemcpyHostToDevice ) );
  }

  void prepareMomenta(fptype energy) override {
    // --- 2a. Fill in momenta of initial state particles on the device
    grambo2toNm0::getMomentaInitial<<<gpublocks, gputhreads>>>( energy, devMomenta.get() );
    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    grambo2toNm0::getMomentaFinal<<<gpublocks, gputhreads>>>( energy, devRnarray.get(), devMomenta.get(), devWeights.get() );
    // --- 2c. CopyDToH Weights
    checkCuda( cudaMemcpy( hstWeights.get(), devWeights.get(),
                           nWeights * sizeof(decltype(hstWeights)::element_type), cudaMemcpyDeviceToHost ) );
    // --- 2d. CopyDToH Momenta
    checkCuda( cudaMemcpy( hstMomenta.get(), devMomenta.get(),
                           nMomenta * sizeof(decltype(hstMomenta)::element_type), cudaMemcpyDeviceToHost ) );
  }

  void runSigmaKin(std::size_t iiter) override {
    // --- 0d. SGoodHel
    if ( iiter == 0 )
    {
      // ... 0d1. Compute good helicity mask on the device
      gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(devMomenta.get(), devMEs.get(), devIsGoodHel.get());
      checkCuda( cudaPeekAtLastError() );
      // ... 0d2. Copy back good helicity mask to the host
      checkCuda( cudaMemcpy( hstIsGoodHel.get(), devIsGoodHel.get(),
                             CPPPROCESS_NCOMB * sizeof(decltype(hstIsGoodHel)::element_type), cudaMemcpyDeviceToHost ) );
      // ... 0d3. Copy back good helicity list to constant memory on the device
      gProc::sigmaKin_setGoodHel(hstIsGoodHel.get());
    }

    // --- 3a. SigmaKin
#ifndef MGONGPU_NSIGHT_DEBUG
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta.get(), devMEs.get());
#else
    gProc::sigmaKin<<<gpublocks, gputhreads, ntpbMAX*sizeof(float)>>>(devMomenta.get(), devMEs.get());
#endif
    checkCuda( cudaPeekAtLastError() );

    // --- 3b. CopyDToH MEs
    checkCuda( cudaMemcpy( hstMEs.get(), devMEs.get(), nMEs * sizeof(decltype(hstMEs)::element_type), cudaMemcpyDeviceToHost ) );
  }

  fptype getMomentum(std::size_t evtNo, unsigned int particle, unsigned int component) const override {
    assert(component < CPPPROCESS_NP4);
    assert(particle  < CPPPROCESS_NPAR);
    const auto page  = evtNo / mgOnGpu::neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % mgOnGpu::neppM; // #event in the current eventpage in this iteration
    return hstMomenta[page * CPPPROCESS_NPAR*CPPPROCESS_NP4*mgOnGpu::neppM +
                      particle * mgOnGpu::neppM*CPPPROCESS_NP4 + component * mgOnGpu::neppM + ieppM];
  };

  fptype getMatrixElement(std::size_t evtNo) const override {
    return hstMEs[evtNo];
  }

};
#endif

// Use two levels of macros to force stringification at the right level
// (see https://gcc.gnu.org/onlinedocs/gcc-3.0.1/cpp_3.html#SEC17 and https://stackoverflow.com/a/3419392)
// Google macro is in https://github.com/google/googletest/blob/master/googletest/include/gtest/gtest-param-test.h
#define TESTID_CPU(s) s##_CPU
#define XTESTID_CPU(s) TESTID_CPU(s)
#define MG_INSTANTIATE_TEST_SUITE_CPU( prefix, test_suite_name )        \
  INSTANTIATE_TEST_SUITE_P( prefix,                                     \
                            test_suite_name,                            \
                            testing::Values( [](){ return new CPUTest; } ) );
#define TESTID_GPU(s) s##_GPU
#define XTESTID_GPU(s) TESTID_GPU(s)
#define MG_INSTANTIATE_TEST_SUITE_GPU( prefix, test_suite_name )        \
  INSTANTIATE_TEST_SUITE_P( prefix,                                     \
                            test_suite_name,                            \
                            testing::Values( [](){ return new CUDATest; } ) );

#if defined MGONGPU_FPTYPE_DOUBLE

#ifdef __CUDACC__
MG_INSTANTIATE_TEST_SUITE_GPU( XTESTID_GPU(MG_EPOCH_PROCESS_ID), MadgraphTestDouble );
#else
MG_INSTANTIATE_TEST_SUITE_CPU( XTESTID_CPU(MG_EPOCH_PROCESS_ID), MadgraphTestDouble );
#endif

#else

#ifdef __CUDACC__
MG_INSTANTIATE_TEST_SUITE_GPU( XTESTID_GPU(MG_EPOCH_PROCESS_ID), MadgraphTestFloat );
#else
MG_INSTANTIATE_TEST_SUITE_CPU( XTESTID_CPU(MG_EPOCH_PROCESS_ID), MadgraphTestFloat );
#endif

#endif

// Add a dummy test just to check the linking (related to issue #143)
/*
#ifdef __CUDACC__
TEST( XTESTID_GPU(MG_EPOCH_PROCESS_ID), dummy ){}
#else
TEST( XTESTID_CPU(MG_EPOCH_PROCESS_ID), dummy ){}
#endif
*/
