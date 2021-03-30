#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#include "MadgraphTest.h"

#include "CommonRandomNumbers.h"
#include "gCPPProcess.h"
#include "Memory.h"
#ifdef __CUDACC__
#include "grambo.cu"
#else
#include "rambo.h"
#endif



struct CUDA_CPU_TestBase : public TestDriverBase<double> {
  static_assert( gputhreads%mgOnGpu::neppR == 0, "ERROR! #threads/block should be a multiple of neppR" );
  static_assert( gputhreads%mgOnGpu::neppM == 0, "ERROR! #threads/block should be a multiple of neppM" );
  static_assert( gputhreads <= mgOnGpu::ntpbMAX, "ERROR! #threads/block should be <= ntpbMAX" );

  const std::size_t nRnarray{ mgOnGpu::np4 * mgOnGpu::nparf * nevt }; // (NB: ASA layout with nevt=npagR*neppR events per iteration)
  const std::size_t nMomenta{ mgOnGpu::np4 * mgOnGpu::npar  * nevt }; // (NB: nevt=npagM*neppM for ASA layouts)
  const std::size_t nWeights{ nevt };
  const std::size_t nMEs    { nevt };

  CUDA_CPU_TestBase() :
  TestDriverBase()
  {
    TestDriverBase::nparticle = mgOnGpu::npar;
  }

};


#ifndef __CUDACC__
struct CPUTest : public CUDA_CPU_TestBase {
  Proc::CPPProcess process;

  // --- 0b. Allocate memory structures
  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  unique_ptr_host<fptype> hstRnarray  { hstMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  unique_ptr_host<fptype> hstMomenta  { hstMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][npar][np4][neppM] (previously was: lp)
  unique_ptr_host<bool  > hstIsGoodHel{ hstMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_host<fptype> hstWeights  { hstMakeUnique<fptype>( nWeights ) };
  unique_ptr_host<fptype> hstMEs      { hstMakeUnique<fptype>( nMEs ) };

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
    std::vector<fptype> rnd = CommonRandomNumbers::generate<fptype>(nRnarray, 1337 + iiter);
    std::copy(rnd.begin(), rnd.end(), hstRnarray.get());
  }


  void prepareMomenta(fptype energy) override {
    // --- 2a. Fill in momenta of initial state particles on the device
    rambo2toNm0::getMomentaInitial( energy, hstMomenta.get(), nevt );

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    rambo2toNm0::getMomentaFinal( energy, hstRnarray.get(), hstMomenta.get(), hstWeights.get(), nevt );
  }


  void runSigmaKin(std::size_t /*iiter*/) override {
    // --- 3a. SigmaKin
    Proc::sigmaKin(hstMomenta.get(), hstMEs.get(), nevt);
  }



  double getMomentum(std::size_t evtNo, unsigned int particle, unsigned int component) const override {
    assert(component < mgOnGpu::np4);
    assert(particle  < mgOnGpu::npar);
    const auto page  = evtNo / mgOnGpu::neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % mgOnGpu::neppM; // #event in the current eventpage in this iteration
    return hstMomenta[page * mgOnGpu::npar*mgOnGpu::np4*mgOnGpu::neppM + particle * mgOnGpu::neppM*mgOnGpu::np4 + component * mgOnGpu::neppM + ieppM];
  };

  double getMatrixElement(std::size_t evtNo) const override {
    return hstMEs[evtNo];
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
  unique_ptr_host<fptype> hstRnarray  { hstMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  unique_ptr_host<fptype> hstMomenta  { hstMakeUnique<fptype>( nMomenta ) }; // AOSOA[npagM][npar][np4][neppM] (previously was: lp)
  unique_ptr_host<bool  > hstIsGoodHel{ hstMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_host<fptype> hstWeights  { hstMakeUnique<fptype>( nWeights ) };
  unique_ptr_host<fptype> hstMEs      { hstMakeUnique<fptype>( nMEs ) };


  unique_ptr_dev<fptype> devRnarray  { devMakeUnique<fptype>( nRnarray ) }; // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  unique_ptr_dev<fptype> devMomenta  { devMakeUnique<fptype>( nMomenta ) }; // (previously was: allMomenta)
  unique_ptr_dev<bool  > devIsGoodHel{ devMakeUnique<bool  >( mgOnGpu::ncomb ) };
  unique_ptr_dev<fptype> devWeights  { devMakeUnique<fptype>( nWeights ) }; // (previously was: meDevPtr)
  unique_ptr_dev<fptype> devMEs      { devMakeUnique<fptype>( nMEs )     }; // (previously was: meDevPtr)

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
    std::vector<fptype> rnd = CommonRandomNumbers::generate<fptype>(nRnarray, 1337 + iiter);
    std::copy(rnd.begin(), rnd.end(), hstRnarray.get());
    checkCuda( cudaMemcpy( devRnarray.get(), hstRnarray.get(), nRnarray * sizeof(decltype(devRnarray)::element_type), cudaMemcpyHostToDevice ) );
  }


  void prepareMomenta(fptype energy) override {
    // --- 2a. Fill in momenta of initial state particles on the device
    grambo2toNm0::getMomentaInitial<<<gpublocks, gputhreads>>>( energy, devMomenta.get() );

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    grambo2toNm0::getMomentaFinal<<<gpublocks, gputhreads>>>( energy, devRnarray.get(), devMomenta.get(), devWeights.get() );

    // --- 2c. CopyDToH Weights
    checkCuda( cudaMemcpy( hstWeights.get(), devWeights.get(), nWeights * sizeof(decltype(hstWeights)::element_type), cudaMemcpyDeviceToHost ) );

    // --- 2d. CopyDToH Momenta
    checkCuda( cudaMemcpy( hstMomenta.get(), devMomenta.get(), nMomenta * sizeof(decltype(hstMomenta)::element_type), cudaMemcpyDeviceToHost ) );
  }


  void runSigmaKin(std::size_t iiter) override {
    // --- 0d. SGoodHel
    if ( iiter == 0 )
    {
      // ... 0d1. Compute good helicity mask on the device
      gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(devMomenta.get(), devIsGoodHel.get());
      checkCuda( cudaPeekAtLastError() );
      // ... 0d2. Copy back good helicity mask to the host
      checkCuda( cudaMemcpy( hstIsGoodHel.get(), devIsGoodHel.get(), mgOnGpu::ncomb * sizeof(decltype(hstIsGoodHel)::element_type), cudaMemcpyDeviceToHost ) );
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


  double getMomentum(std::size_t evtNo, unsigned int particle, unsigned int component) const override {
    assert(component < mgOnGpu::np4);
    assert(particle  < mgOnGpu::npar);
    const auto page  = evtNo / mgOnGpu::neppM; // #eventpage in this iteration
    const auto ieppM = evtNo % mgOnGpu::neppM; // #event in the current eventpage in this iteration
    return hstMomenta[page * mgOnGpu::npar*mgOnGpu::np4*mgOnGpu::neppM + particle * mgOnGpu::neppM*mgOnGpu::np4 + component * mgOnGpu::neppM + ieppM];
  };

  double getMatrixElement(std::size_t evtNo) const override {
    return hstMEs[evtNo];
  }
};
#endif


#ifdef __CUDACC__
INSTANTIATE_TEST_SUITE_P(EP2_CUDA_GPU, MadgraphTestDouble,
    testing::Values( [](){ return new CUDATest; } )
);
#else
INSTANTIATE_TEST_SUITE_P(EP2_CUDA_CPU, MadgraphTestDouble,
    testing::Values([](){ return new CPUTest; })
);
#endif


