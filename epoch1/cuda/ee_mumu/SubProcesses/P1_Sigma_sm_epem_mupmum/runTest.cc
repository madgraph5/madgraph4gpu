#include "epoch_process_id.h"

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "MadgraphTest.h"

#include "CommonRandomNumbers.h"
#include "CPPProcess.h"
#include "Memory.h"
#ifdef __CUDACC__
#include "grambo.cu"
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
using unique_ptr_host = std::unique_ptr<T[], CppHstDeleter<T>>;
#endif

struct CUDA_CPU_TestBase : public TestDriverBase {

  static_assert( gputhreads%mgOnGpu::neppR == 0, "ERROR! #threads/block should be a multiple of neppR" );
  static_assert( gputhreads%mgOnGpu::neppM == 0, "ERROR! #threads/block should be a multiple of neppM" );
  static_assert( gputhreads <= mgOnGpu::ntpbMAX, "ERROR! #threads/block should be <= ntpbMAX" );

  const std::size_t nRnarray{ mgOnGpu::np4 * mgOnGpu::nparf * nevt }; // AOSOA layout with nevt=npagR*neppR events per iteration
  const std::size_t nMomenta{ mgOnGpu::np4 * mgOnGpu::npar  * nevt }; // AOSOA layout with nevt=npagM*neppM events per iteration
  const std::size_t nWeights{ nevt };
  const std::size_t nMEs    { nevt };

  CUDA_CPU_TestBase() :
  TestDriverBase()
  {  }

};

#ifndef __CUDACC__
struct CPUTest : public CUDA_CPU_TestBase {

  Proc::CPPProcess process;

  // --- 0b. Allocate memory structures
  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  unique_ptr_host<fptype> hstRnarray  { hstMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR]
  unique_ptr_host<fptype> hstMomenta  { hstMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][npar][np4][neppM]
  unique_ptr_host<bool  > hstIsGoodHel{ hstMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_host<fptype> hstWeights  { hstMakeUnique<fptype>( nWeights ) };
  unique_ptr_host<fptype> hstMEs      { hstMakeUnique<fptype>( nMEs ) }; // ARRAY[nevt]

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
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    using mgOnGpu::neppM;
    assert(component < np4);
    assert(particle  < npar);
    const auto ipagM = evtNo / neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % neppM; // #event in the current eventpage in this iteration
    return hstMomenta[ipagM*npar*np4*neppM + particle*np4*neppM + component*neppM + ieppM];
    /*
#ifndef MGONGPU_CPPSIMD
    return hstMomenta[ipagM*npar*np4*neppM + particle*np4*neppM + component*neppM + ieppM];
#else
    return hstMomenta[ipagM*npar*np4 + particle*np4 + component][ieppM];
#endif
    */
  };

  fptype getMatrixElement(std::size_t ievt) const override {
    return hstMEs[ievt];
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
  unique_ptr_host<fptype> hstRnarray  { hstMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR] (nevt=npagR*neppR)
  unique_ptr_host<fptype> hstMomenta  { hstMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][npar][np4][neppM] (nevt=npagM*neppM)
  unique_ptr_host<bool  > hstIsGoodHel{ hstMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_host<fptype> hstWeights  { hstMakeUnique<fptype>( nWeights ) };
  unique_ptr_host<fptype> hstMEs      { hstMakeUnique<fptype>( nMEs ) }; // ARRAY[nevt]

  unique_ptr_dev<fptype> devRnarray  { devMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR] (nevt=npagR*neppR)
  unique_ptr_dev<fptype> devMomenta  { devMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][npar][np4][neppM] (nevt=npagM*neppM)
  unique_ptr_dev<bool  > devIsGoodHel{ devMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_dev<fptype> devWeights  { devMakeUnique<fptype>( nWeights ) };
  unique_ptr_dev<fptype> devMEs      { devMakeUnique<fptype>( nMEs ) }; // ARRAY[nevt]

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
                             mgOnGpu::ncomb * sizeof(decltype(hstIsGoodHel)::element_type), cudaMemcpyDeviceToHost ) );
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
    assert(component < mgOnGpu::np4);
    assert(particle  < mgOnGpu::npar);
    const auto page  = evtNo / mgOnGpu::neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % mgOnGpu::neppM; // #event in the current eventpage in this iteration
    return hstMomenta[page * mgOnGpu::npar*mgOnGpu::np4*mgOnGpu::neppM +
                      particle * mgOnGpu::neppM*mgOnGpu::np4 + component * mgOnGpu::neppM + ieppM];
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
                            testing::Values( new CPUTest ) );
#define TESTID_GPU(s) s##_GPU
#define XTESTID_GPU(s) TESTID_GPU(s)
#define MG_INSTANTIATE_TEST_SUITE_GPU( prefix, test_suite_name )        \
  INSTANTIATE_TEST_SUITE_P( prefix,                                     \
                            test_suite_name,                            \
                            testing::Values( new CUDATest ) );

#ifdef __CUDACC__
MG_INSTANTIATE_TEST_SUITE_GPU( XTESTID_GPU(MG_EPOCH_PROCESS_ID), MadgraphTest );
#else
MG_INSTANTIATE_TEST_SUITE_CPU( XTESTID_CPU(MG_EPOCH_PROCESS_ID), MadgraphTest );
#endif
