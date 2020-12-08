#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <unistd.h>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#ifdef __CUDACC__
#include "grambo.cu"
#else
#include "rambo.h"
#endif

#ifdef MGONGPU_COMMONRAND_ONHOST
#include "CommonRandomNumbers.h"
#endif
#include "CPPProcess.h"
#include "timermap.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return (int)strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only controls if nan's emit warnings" << std::endl;
  return ret;
}

#ifdef __CUDACC__
template<typename T = fptype>
struct CudaDevDeleter {
  void operator()(T* mem) {
    checkCuda( cudaFree( mem ) );
  }
};
template<typename T = fptype>
std::unique_ptr<T, CudaDevDeleter<T>> devMakeUnique(std::size_t N) {
  T* tmp = nullptr;
  checkCuda( cudaMalloc( &tmp, N * sizeof(T) ) );
  return std::unique_ptr<T, CudaDevDeleter<T>>{ tmp };
}
template<typename T = fptype>
struct CudaHstDeleter {
  void operator()(T* mem) {
    checkCuda( cudaFreeHost( mem ) );
  }
};
template<typename T = fptype>
std::unique_ptr<T[], CudaHstDeleter<T>> hstMakeUnique(std::size_t N) {
  T* tmp = nullptr;
  checkCuda( cudaMallocHost( &tmp, N * sizeof(T) ) );
  return std::unique_ptr<T[], CudaHstDeleter<T>>{ tmp };
};
#else
template<typename T = fptype>
std::unique_ptr<T[]> hstMakeUnique(std::size_t N) { return std::unique_ptr<T[]>{ new T[N]() }; };
#endif

int main(int argc, char **argv)
{
  sycl::queue q_ct1{ sycl::cpu_selector{} };
  // READ COMMAND LINE ARGUMENTS
  bool verbose = false;
  bool debug = false;
  bool perf = false;
  bool json = false;
  int niter = 0;
  int gpublocks = 1;
  int gputhreads = 32;
  int jsondate = 0;
  int jsonrun = 0;
  int numvec[5] = {0,0,0,0,0};
  int nnum = 0;

  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
    else if (strcmp(argv[argn], "--performance") == 0 ||
             strcmp(argv[argn], "-p") == 0)
      perf = true;
    else if (strcmp(argv[argn], "--json") == 0 ||
             strcmp(argv[argn], "-j") == 0)
      json = true;
    else if (is_number(argv[argn]) && nnum<5)
      numvec[nnum++] = atoi(argv[argn]);
    else
      return usage(argv[0]);
  }

  if (nnum == 3 || nnum == 5) {
    gpublocks = numvec[0];
    gputhreads = numvec[1];
    niter = numvec[2];
    if (nnum == 5){
      jsondate = numvec[3];
      jsonrun = numvec[4];
    }
  } else if (nnum == 1) {
    niter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (niter == 0)
    return usage(argv[0]);

  const int neppR = mgOnGpu::neppR; // ASA layout: constant at compile-time
  if ( gputhreads%neppR != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of neppR=" << neppR << std::endl;
    return usage(argv[0]);
  }

  const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
  if ( gputhreads%neppM != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of neppM=" << neppM << std::endl;
    return usage(argv[0]);
  }

  using mgOnGpu::ntpbMAX;
  if ( gputhreads > ntpbMAX )
  {
    std::cout << "ERROR! #threads/block should be <= " << ntpbMAX << std::endl;
    return usage(argv[0]);
  }

  const int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  const int nevt = ndim; // number of events in one iteration == number of GPU threads
  const int nevtALL = niter*nevt; // total number of ALL events in all iterations

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // *** START THE NEW TIMERS ***
  mgOnGpu::TimerMap timermap;

  // === STEP 0 - INITIALISE

#ifdef __CUDACC__
  // --- 00. Initialise cuda (call cudaFree to ease cuda profile analysis)
  const std::string cdfrKey = "00 CudaFree";
  timermap.start( cdfrKey );
  //std::cout << "Calling cudaFree... " << std::endl;
  checkCuda( cudaFree( 0 ) ); // SLOW!
  //std::cout << "Calling cudaFree... done" << std::endl;

  // --- Book the tear down at the end of main:
  struct CudaTearDown {
    CudaTearDown(bool print) : _print(print) { }
    ~CudaTearDown() {
      if ( _print ) std::cout << "Calling cudaDeviceReset()." << std::endl;
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
    bool _print{false};
  } cudaTearDown(debug);
#endif

  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
#ifdef __CUDACC__
  gProc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#else
  Proc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#endif

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");
  const fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
  //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
  const int meGeVexponent = -(2 * process.nexternal - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  using mgOnGpu::np4;
  using mgOnGpu::nparf;
  using mgOnGpu::npar;
  using mgOnGpu::ncomb; // Number of helicity combinations
  const unsigned int nRnarray = np4*nparf*nevt; // (NB: ASA layout with nevt=npagR*neppR events per iteration)
  const unsigned int nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for ASA layouts)
  const unsigned int nWeights = nevt;
  const unsigned int nMEs     = nevt;

#if defined MGONGPU_CURAND_ONHOST or defined MGONGPU_COMMONRAND_ONHOST or not defined __CUDACC__
  auto hstRnarray   = hstMakeUnique<fptype>( nRnarray ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
#endif
  sycl::buffer<fptype, 1, sycl::buffer_allocator> Momenta_buffer{sycl::range<1>{nMomenta}};
  sycl::buffer<bool, 1, sycl::buffer_allocator> IsGoodHel_buffer{sycl::range<1>{ncomb}};
  sycl::buffer<fptype, 1, sycl::buffer_allocator> Weights_buffer{sycl::range<1>{nWeights}};
  sycl::buffer<fptype, 1, sycl::buffer_allocator> MEs_buffer{sycl::range<1>{nMEs}};
  sycl::buffer<int, 2, sycl::buffer_allocator> cHel_buffer{sycl::range<2>{ncomb,npar}};
  sycl::buffer<int, 1, sycl::buffer_allocator> cNGoodHel_buffer{sycl::range<1>{1}};
  sycl::buffer<int, 1, sycl::buffer_allocator> cGoodHel_buffer{sycl::range<1>{ncomb}};
#ifdef __CUDACC__
  auto devRnarray   = devMakeUnique<fptype>( nRnarray ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  auto devMomenta   = devMakeUnique<fptype>( nMomenta ); // (previously was: allMomenta)
  auto devIsGoodHel = devMakeUnique<bool  >( ncomb );
  auto devWeights   = devMakeUnique<fptype>( nWeights ); // (previously was: meDevPtr)
  auto devMEs       = devMakeUnique<fptype>( nMEs ); // (previously was: meDevPtr)

#if defined MGONGPU_CURAND_ONHOST or defined MGONGPU_COMMONRAND_ONHOST
  const int nbytesRnarray = nRnarray * sizeof(fptype);
#endif
  const int nbytesMomenta = nMomenta * sizeof(fptype);
  const int nbytesIsGoodHel = ncomb * sizeof(bool);
  const int nbytesWeights = nWeights * sizeof(fptype);
  const int nbytesMEs = nMEs * sizeof(fptype);
#endif

  std::unique_ptr<double[]> genrtimes( new double[niter] );
  std::unique_ptr<double[]> rambtimes( new double[niter] );
  std::unique_ptr<double[]> wavetimes( new double[niter] );
  std::unique_ptr<fptype[]> matrixelementALL( new fptype[nevtALL] ); // FIXME: assume process.nprocesses == 1
  std::unique_ptr<fptype[]> weightALL( new fptype[nevtALL] );

  // --- 0c. Create curand or common generator
  const std::string cgenKey = "0c GenCreat";
  timermap.start( cgenKey );
#ifdef MGONGPU_COMMONRAND_ONHOST
  std::vector<std::promise<std::vector<fptype>>> commonRandomPromises;
  CommonRandomNumbers::startGenerateAsync(commonRandomPromises, nRnarray, niter);
#else
  curandGenerator_t rnGen;
#ifdef __CUDACC__
  grambo2toNm0::createGenerator( &rnGen );
#else
  rambo2toNm0::createGenerator( &rnGen );
#endif
#endif

  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************

  for (int iiter = 0; iiter < niter; ++iiter)
  {
    //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;

    // === STEP 1 OF 3

    // *** START THE OLD-STYLE TIMER FOR RANDOM GEN ***
    double genrtime = 0;

#if defined MGONGPU_CURAND_ONHOST or defined MGONGPU_CURAND_ONDEVICE
    // --- 1a. Seed curand generator (to get same results on host and device)
    // [NB This should not be necessary using the host API: "Generation functions
    // can be called multiple times on the same generator to generate successive
    // blocks of results. For pseudorandom generators, multiple calls to generation
    // functions will yield the same result as a single call with a large size."]
    const std::string sgenKey = "1a GenSeed ";
    timermap.start( sgenKey );
    const unsigned long long seed = 20200805;
#ifdef __CUDACC__
    grambo2toNm0::seedGenerator( rnGen, seed+iiter );
#else
    rambo2toNm0::seedGenerator( rnGen, seed+iiter );
#endif
    genrtime += timermap.stop();
#endif

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );
#ifdef MGONGPU_COMMONRAND_ONHOST
    std::vector<fptype> commonRnd = commonRandomPromises[iiter].get_future().get();
    assert( nRnarray == static_cast<int>( commonRnd.size() ) );
    // NB (PR #45): memcpy is strictly needed only in CUDA (copy to pinned memory), but keep it also in C++ for consistency
    memcpy( hstRnarray.get(), commonRnd.data(), nRnarray * sizeof(fptype) );
#elif defined __CUDACC__
#ifdef MGONGPU_CURAND_ONDEVICE
    grambo2toNm0::generateRnarray( rnGen, devRnarray.get(), nevt );
#elif defined MGONGPU_CURAND_ONHOST
    grambo2toNm0::generateRnarray( rnGen, hstRnarray.get(), nevt );
#endif
#else
    rambo2toNm0::generateRnarray( rnGen, hstRnarray.get(), nevt );
#endif
    //std::cout << "Got random numbers" << std::endl;

#ifdef __CUDACC__
#ifndef MGONGPU_CURAND_ONDEVICE
    // --- 1c. Copy rnarray from host to device
    const std::string htodKey = "1c CpHTDrnd";
    genrtime += timermap.start( htodKey );
    // NB (PR #45): this cudaMemcpy would involve an intermediate memcpy to pinned memory, if hstRnarray was not already cudaMalloc'ed
    checkCuda( cudaMemcpy( devRnarray.get(), hstRnarray.get(), nbytesRnarray, cudaMemcpyHostToDevice ) );
#endif
#endif

    // *** STOP THE OLD-STYLE TIMER FOR RANDOM GEN ***
    genrtime += timermap.stop();

    // === STEP 2 OF 3
    // Fill in particle momenta for each of nevt events on the device

    // *** START THE OLD-STYLE TIMER FOR RAMBO ***
    double rambtime = 0;

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );
#ifdef SYCL_LANGUAGE_VERSION
    {
      q_ct1.submit([&](sycl::handler &cgh) {
        auto devMomenta_acc = Momenta_buffer.get_access<sycl::access::mode::read_write>(cgh);
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
              auto devMomenta_ptr = devMomenta_acc.get_pointer();
              rambo2toNm0::getMomentaInitial(energy, devMomenta_ptr, item_ct1);
            });
      });
      q_ct1.wait();
    }
#else
    rambo2toNm0::getMomentaInitial( energy, hstMomenta.get(), nevt );
#endif
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );
#ifdef SYCL_LANGUAGE_VERSION
    {
      sycl::buffer Rnarray_buffer{hstRnarray.get(), sycl::range{nRnarray}};
      q_ct1.submit([&](sycl::handler &cgh) {
        auto devRnarray_acc = Rnarray_buffer.get_access<sycl::access::mode::read>(cgh);
        auto devMomenta_acc = Momenta_buffer.get_access<sycl::access::mode::read_write>(cgh);
        auto devWeights_acc = Weights_buffer.get_access<sycl::access::mode::write>(cgh);
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
              auto devRnarray_ptr = devRnarray_acc.get_pointer();
              auto devWeights_ptr = devWeights_acc.get_pointer();
              auto devMomenta_ptr = devMomenta_acc.get_pointer();
              rambo2toNm0::getMomentaFinal(energy, devRnarray_ptr,
                                            devMomenta_ptr,
                                            devWeights_ptr, item_ct1);
            });
      });
    q_ct1.wait();

    }
#else
    rambo2toNm0::getMomentaFinal( energy, hstRnarray.get(), hstMomenta.get(), hstWeights.get(), nevt );
#endif
    //std::cout << "Got final momenta" << std::endl;

#ifdef __CUDACC__
    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    rambtime += timermap.start( cwgtKey );
    checkCuda( cudaMemcpy( hstWeights.get(), devWeights.get(), nbytesWeights, cudaMemcpyDeviceToHost ) );


    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    rambtime += timermap.start( cmomKey );
    checkCuda( cudaMemcpy( hstMomenta.get(), devMomenta.get(), nbytesMomenta, cudaMemcpyDeviceToHost ) );

#endif

    // *** STOP THE OLD-STYLE TIMER FOR RAMBO ***
    rambtime += timermap.stop();

    // === STEP 3 OF 3
    // Evaluate matrix elements for all nevt events
    // 0d. (Only on the first iteration) Get good helicities [renamed as 0d: this is initialisation!]
    // 3a. Evaluate MEs on the device
    // 3b. Copy MEs back from device to host

    // --- 0d. SGoodHel
#ifdef SYCL_LANGUAGE_VERSION
    if ( iiter == 0 )
    {
      const std::string ghelKey = "0d SGoodHel";
      timermap.start( ghelKey );
      // ... 0d1. Compute good helicity mask on the device
      {
        q_ct1.submit([&](sycl::handler &cgh) {
          auto devcHel_acc = cHel_buffer.get_access<sycl::access::mode::read_write>(cgh);
          auto devMomenta_acc = Momenta_buffer.get_access<sycl::access::mode::read_write>(cgh);
          auto devIsGoodHel_acc = IsGoodHel_buffer.get_access<sycl::access::mode::read_write>(cgh);
          cgh.parallel_for(
              sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                    sycl::range<3>(1, 1, gputhreads),
                                sycl::range<3>(1, 1, gputhreads)),
              [=](sycl::nd_item<3> item_ct1) {
                auto devMomenta_ptr = devMomenta_acc.get_pointer();
                auto devIsGoodHel_ptr = devIsGoodHel_acc.get_pointer();
                Proc::sigmaKin_getGoodHel(
                devMomenta_ptr, devIsGoodHel_ptr, item_ct1, devcHel_acc);
              });
        });
      q_ct1.wait();
      }
      // ... 0d2. Copy back good helicity mask to the host
      sycl::host_accessor IsGoodHel_acc{IsGoodHel_buffer};
      sycl::host_accessor cNGoodHel_acc{cNGoodHel_buffer};
      sycl::host_accessor cGoodHel_acc{cGoodHel_buffer};
      auto IsGoodHel_ptr = IsGoodHel_acc.get_pointer();
      auto cNGoodHel_ptr = cNGoodHel_acc.get_pointer();
      auto cGoodHel_ptr = cGoodHel_acc.get_pointer();
      Proc::sigmaKin_setGoodHel(IsGoodHel_ptr, cNGoodHel_ptr, cGoodHel_ptr) ;
      // ... 0d3. Copy back good helicity list to constant memory on the device
    }
#endif

    // *** START THE OLD TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    double wavetime = 0;

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
#ifdef SYCL_LANGUAGE_VERSION
#ifndef MGONGPU_NSIGHT_DEBUG
    {
      q_ct1.submit([&](sycl::handler &cgh) {
        auto cHel_acc = cHel_buffer.get_access<sycl::access::mode::read_write>(cgh);
        auto cNGoodHel_acc = cNGoodHel_buffer.get_access<sycl::access::mode::read_write>(cgh);
        auto cGoodHel_acc = cGoodHel_buffer.get_access<sycl::access::mode::read_write>(cgh);
        auto devMomenta_acc = Momenta_buffer.get_access<sycl::access::mode::read_write>(cgh);
        auto devMEs_acc = MEs_buffer.get_access<sycl::access::mode::read_write>(cgh);
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
                auto devMomenta_ptr = devMomenta_acc.get_pointer();
                auto devMEs_ptr = devMEs_acc.get_pointer();
                auto cNGoodHel_ptr = cNGoodHel_acc.get_pointer();
                auto cGoodHel_ptr = cGoodHel_acc.get_pointer();
                Proc::sigmaKin(
                devMomenta_ptr, devMEs_ptr, item_ct1, cHel_acc, cNGoodHel_ptr, cGoodHel_ptr);
              });
      });
    q_ct1.wait();
    }
#else
    gProc::sigmaKin<<<gpublocks, gputhreads, ntpbMAX*sizeof(float)>>>(devMomenta.get(), devMEs.get());
#endif
#else
    Proc::sigmaKin(hstMomenta.get(), hstMEs.get(), nevt);
#endif

#ifdef SYCL_LANGUAGE_VERSION
    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    wavetime += timermap.start( cmesKey );
#endif

    // *** STOP THE OLD TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wavetime += timermap.stop();

    // === STEP 4 FINALISE LOOP
    // --- 4a Dump within the loop
    const std::string loopKey = "4a DumpLoop";
    timermap.start(loopKey);
    genrtimes[iiter] = genrtime;
    rambtimes[iiter] = rambtime;
    wavetimes[iiter] = wavetime;

    if (verbose)
    {
      std::cout << "***********************************************************************" << std::endl
                << "Iteration #" << iiter+1 << " of " << niter << std::endl;
      if (perf) std::cout << "Wave function time: " << wavetime << std::endl;
    }

    sycl::host_accessor Weights{Weights_buffer};
    sycl::host_accessor Momenta{Momenta_buffer};
    sycl::host_accessor MEs{MEs_buffer};
    for (int ievt = 0; ievt < nevt; ++ievt) // Loop over all events in this iteration
    {
      if (verbose)
      {
        // Display momenta
        const int ipagM = ievt/neppM; // #eventpage in this iteration
        const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
        std::cout << "Momenta:" << std::endl;
        for (int ipar = 0; ipar < npar; ipar++)
        {
          // NB: 'setw' affects only the next field (of any type)
          std::cout << std::scientific // fixed format: affects all floats (default precision: 6)
                    << std::setw(4) << ipar + 1
                    << std::setw(14) << Momenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 0*neppM + ieppM] // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << Momenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 1*neppM + ieppM] // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << Momenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 2*neppM + ieppM] // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << Momenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 3*neppM + ieppM] // AOSOA[ipagM][ipar][3][ieppM]
                    << std::endl
                    << std::defaultfloat; // default format: affects all floats
        }
        std::cout << std::string(80, '-') << std::endl;
        // Display matrix elements
        std::cout << " Matrix element = "
                  << MEs[ievt] << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(80, '-') << std::endl;
      }
      // Fill the arrays with ALL MEs and weights
      matrixelementALL[iiter*nevt + ievt] = MEs[ievt]; // FIXME: assume process.nprocesses == 1
      weightALL[iiter*nevt + ievt] = Weights[ievt];
    }

    if (!(verbose || debug || perf))
    {
      std::cout << ".";
    }

  }

  // **************************************
  // *** END MAIN LOOP ON #ITERATIONS ***
  // **************************************

  // === STEP 8 ANALYSIS
  // --- 8a Analysis: compute stats after the loop
  const std::string statKey = "8a CompStat";
  timermap.start(statKey);

  double sumgtim = 0;
  double sqsgtim = 0;
  double mingtim = genrtimes[0];
  double maxgtim = genrtimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumgtim += genrtimes[iiter];
    sqsgtim += genrtimes[iiter]*genrtimes[iiter];
    mingtim = std::min( mingtim, genrtimes[iiter] );
    maxgtim = std::max( maxgtim, genrtimes[iiter] );
  }

  double sumrtim = 0;
  double sqsrtim = 0;
  double minrtim = rambtimes[0];
  double maxrtim = rambtimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumrtim += rambtimes[iiter];
    sqsrtim += rambtimes[iiter]*rambtimes[iiter];
    minrtim = std::min( minrtim, rambtimes[iiter] );
    maxrtim = std::max( maxrtim, rambtimes[iiter] );
  }

  double sumwtim = 0;
  double sqswtim = 0;
  double minwtim = wavetimes[0];
  double maxwtim = wavetimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumwtim += wavetimes[iiter];
    sqswtim += wavetimes[iiter]*wavetimes[iiter];
    minwtim = std::min( minwtim, wavetimes[iiter] );
    maxwtim = std::max( maxwtim, wavetimes[iiter] );
  }
  double meanwtim = sumwtim / niter;
  //double stdwtim = std::sqrt( sqswtim / niter - meanwtim * meanwtim );

  int nnan = 0;
  double minelem = matrixelementALL[0];
  double maxelem = matrixelementALL[0];
  double minweig = weightALL[0];
  double maxweig = weightALL[0];
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute min/max
    if ( std::isnan( matrixelementALL[ievtALL] ) )
    {
      if ( debug ) // only printed out with "-p -d" (matrixelementALL is not filled without -p)
        std::cout << "WARNING! ME[" << ievtALL << "} is nan" << std::endl;
      nnan++;
      continue;
    }
    minelem = std::min( minelem, (double)matrixelementALL[ievtALL] );
    maxelem = std::max( maxelem, (double)matrixelementALL[ievtALL] );
    minweig = std::min( minweig, (double)weightALL[ievtALL] );
    maxweig = std::max( maxweig, (double)weightALL[ievtALL] );
  }
  double sumelemdiff = 0;
  double sumweigdiff = 0;
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute mean from the sum of diff to min
    if ( std::isnan( matrixelementALL[ievtALL] ) ) continue;
    sumelemdiff += ( matrixelementALL[ievtALL] - minelem );
    sumweigdiff += ( weightALL[ievtALL] - minweig );
  }
  double meanelem = minelem + sumelemdiff / ( nevtALL - nnan );
  double meanweig = minweig + sumweigdiff / ( nevtALL - nnan );
  double sqselemdiff = 0;
  double sqsweigdiff = 0;
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute stddev from the squared sum of diff to mean
    if ( std::isnan( matrixelementALL[ievtALL] ) ) continue;
    sqselemdiff += std::pow( matrixelementALL[ievtALL] - meanelem, 2 );
    sqsweigdiff += std::pow( weightALL[ievtALL] - meanweig, 2 );
  }
  double stdelem = std::sqrt( sqselemdiff / ( nevtALL - nnan ) );
  double stdweig = std::sqrt( sqsweigdiff / ( nevtALL - nnan ) );

  // === STEP 9 FINALISE
  // --- 9a. Destroy curand generator
  const std::string dgenKey = "9a GenDestr";
  timermap.start( dgenKey );
#ifndef MGONGPU_COMMONRAND_ONHOST
#ifdef __CUDACC__
  grambo2toNm0::destroyGenerator( rnGen );
#else
  rambo2toNm0::destroyGenerator( rnGen );
#endif
#endif

  // --- 9b Dump to screen
  const std::string dumpKey = "9b DumpScrn";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
    std::cout << "***********************************************************************" << std::endl
              << "NumBlocksPerGrid           = " << gpublocks << std::endl
              << "NumThreadsPerBlock         = " << gputhreads << std::endl
              << "NumIterations              = " << niter << std::endl
              << "-----------------------------------------------------------------------" << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
              << "FP precision               = DOUBLE (nan=" << nnan << ")" << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
              << "FP precision               = FLOAT (nan=" << nnan << ")" << std::endl
#endif
#ifdef SYCL_LANGUAGE_VERSION
#if defined MGONGPU_CXTYPE_CUCOMPLEX
              << "Complex type               = CUCOMPLEX" << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
              << "Complex type               = THRUST::COMPLEX" << std::endl
#endif
#else
              << "Complex type               = STD::COMPLEX" << std::endl
#endif
              << "RanNumb memory layout      = AOSOA[" << neppR << "]"
              << ( neppR == 1 ? " == AOS" : "" ) << std::endl
              << "Momenta memory layout      = AOSOA[" << neppM << "]"
              << ( neppM == 1 ? " == AOS" : "" ) << std::endl
#ifdef SYCL_LANGUAGE_VERSION
              << "Wavefunction GPU memory    = LOCAL" << std::endl
#endif
#ifdef SYCL_LANGUAGE_VERSION
#if defined MGONGPU_COMMONRAND_ONHOST
              << "Random number generation   = COMMON RANDOM HOST (CUDA code)" << std::endl
#elif defined MGONGPU_CURAND_ONDEVICE
              << "Random number generation   = CURAND DEVICE (CUDA code)" << std::endl
#elif defined MGONGPU_CURAND_ONHOST
              << "Random number generation   = CURAND HOST (CUDA code)" << std::endl
#endif
#else
#if defined MGONGPU_COMMONRAND_ONHOST
              << "Random number generation   = COMMON RANDOM (C++ code)" << std::endl
#else
              << "Random number generation   = CURAND (C++ code)" << std::endl
#endif
#endif
              << "-----------------------------------------------------------------------" << std::endl
              << "NumberOfEntries            = " << niter << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "TotalTime[Rnd+Rmb+ME] (123)= ( " << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo+ME]    (23)= ( " << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[RndNumGen]    (1)= ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo]        (2)= ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[MatrixElems]  (3)= ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "MeanTimeInMatrixElems      = ( " << meanwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "[Min,Max]TimeInMatrixElems = [ " << minwtim
              << " ,  " << maxwtim << " ]  sec" << std::endl
      //<< "StdDevTimeInWaveFuncs      = ( " << stdwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "-----------------------------------------------------------------------" << std::endl
      //<< "ProcessID:                 = " << getpid() << std::endl
      //<< "NProcesses                 = " << process.nprocesses << std::endl
              << "TotalEventsComputed        = " << nevtALL << std::endl
              << "EvtsPerSec[Rnd+Rmb+ME](123)= ( " << nevtALL/(sumgtim+sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[Rmb+ME]     (23)= ( " << nevtALL/(sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[RndNumbGen]   (1)= ( " << nevtALL/sumgtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[Rambo]        (2)= ( " << nevtALL/sumrtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[MatrixElems] (3)= ( " << nevtALL/sumwtim
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << std::defaultfloat; // default format: affects all floats
    std::cout << "***********************************************************************" << std::endl
              << "NumMatrixElements(notNan)  = " << nevtALL - nnan << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "MeanMatrixElemValue        = ( " << meanelem
              << " +- " << stdelem/sqrt(nevtALL - nnan) << " )  GeV^" << meGeVexponent << std::endl // standard error
              << "[Min,Max]MatrixElemValue   = [ " << minelem
              << " ,  " << maxelem << " ]  GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue      = ( " << stdelem << std::string(16, ' ') << " )  GeV^" << meGeVexponent << std::endl
              << "MeanWeight                 = ( " << meanweig
              << " +- " << stdweig/sqrt(nevtALL - nnan) << " )" << std::endl // standard error
              << "[Min,Max]Weight            = [ " << minweig
              << " ,  " << maxweig << " ]" << std::endl
              << "StdDevWeight               = ( " << stdweig << std::string(16, ' ') << " )" << std::endl
              << std::defaultfloat; // default format: affects all floats
  }

  // --- 9c Dump to json
  const std::string jsonKey = "9c DumpJson";
  timermap.start(jsonKey);

  if(json)
  {
    std::string jsonFileName = std::to_string(jsondate) + "-perf-test-run" + std::to_string(jsonrun) + ".json";
    jsonFileName = "./perf/data/" + jsonFileName;

    //Checks if file exists
    std::ifstream fileCheck;
    bool fileExists = false;
    fileCheck.open(jsonFileName);
    if(fileCheck){
      fileExists = true;
      fileCheck.close();
    }

    std::ofstream jsonFile;
    jsonFile.open(jsonFileName, std::ios_base::app);
    if(!fileExists){
      jsonFile << "[" << std::endl;
    }
    else{
      //deleting the last bracket and outputting a ", "
      std::string temp = "truncate -s-1 " + jsonFileName;
      const char *command = temp.c_str();
      system(command);
      jsonFile << ", " << std::endl;
    }

    jsonFile << "{" << std::endl
             << "\"NumIterations\": " << niter << ", " << std::endl
             << "\"NumThreadsPerBlock\": " << gputhreads << ", " << std::endl
             << "\"NumBlocksPerGrid\": " << gpublocks << ", " << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
             << "\"FP precision\": "
             << "\"DOUBLE (nan=" << nnan << ")\"," << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
             << "\"FP precision\": " << "FLOAT (nan=" << nnan << ")," << std::endl
#endif
             << "\"Complex type\": "
#ifdef __CUDACC__
#if defined MGONGPU_CXTYPE_CUCOMPLEX
             << "\"CUCOMPLEX\"," << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
             << "\"THRUST::COMPLEX\"," << std::endl
#endif
#else
             << "\"STD::COMPLEX\"," << std::endl
#endif
             << "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\""
             << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Momenta memory layout\": " << "\"AOSOA[" << neppM << "]\""
             << ( neppM == 1 ? " == AOS" : "" ) << ", " << std::endl
#ifdef __CUDACC__
             << "\"Wavefunction GPU memory\": " << "\"LOCAL\"," << std::endl
#endif
             << "\"Curand generation\": "
#ifdef __CUDACC__
#if defined MGONGPU_COMMONRAND_ONHOST
             << "\"COMMON RANDOM HOST (CUDA code)\"," << std::endl
#elif defined MGONGPU_CURAND_ONDEVICE
             << "\"CURAND DEVICE (CUDA code)\"," << std::endl
#elif defined MGONGPU_CURAND_ONHOST
             << "\"CURAND HOST (CUDA code)\"," << std::endl
#endif
#else
#if defined MGONGPU_COMMONRAND_ONHOST
             << "\"COMMON RANDOM (C++ code)\"," << std::endl
#else
             << "\"CURAND (C++ code)\"," << std::endl
#endif
#endif
             << "\"NumberOfEntries\": " << niter << "," << std::endl
      //<< std::scientific // Not sure about this

             << "\"TotalTime[Rnd+Rmb+ME] (123)\": \""
             << std::to_string(sumgtim+sumrtim+sumwtim) << " sec\","
             << std::endl
             << "\"TotalTime[Rambo+ME] (23)\": \""
             << std::to_string(sumrtim+sumwtim) << " sec\"," << std::endl
             << "\"TotalTime[RndNumGen] (1)\": \""
             << std::to_string(sumgtim) << " sec\"," << std::endl
             << "\"TotalTime[Rambo] (2)\": \""
             << std::to_string(sumrtim) << " sec\"," << std::endl
             << "\"TotalTime[MatrixElems] (3)\": \""
             << std::to_string(sumwtim) << " sec\"," << std::endl
             << "\"MeanTimeInMatrixElems\": \""
             << std::to_string(meanwtim) << " sec\"," << std::endl
             << "\"MinTimeInMatrixElems\": \""
             << std::to_string(minwtim) << " sec\"," << std::endl
             << "\"MaxTimeInMatrixElems\": \""
             << std::to_string(maxwtim) << " sec\"," << std::endl
      //<< "ProcessID:                = " << getpid() << std::endl
      //<< "NProcesses                = " << process.nprocesses << std::endl
             << "\"TotalEventsComputed\": " << nevtALL << "," << std::endl

             << "\"EvtsPerSec[Rnd+Rmb+ME](123)\": \""
             << std::to_string(nevtALL/(sumgtim+sumrtim+sumwtim))
             << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[Rmb+ME] (23)\": \""
             << std::to_string(nevtALL/(sumrtim+sumwtim)) << " sec^-1\","
             << std::endl
             << "\"EvtsPerSec[MatrixElems] (3)\": \""
             << std::to_string(nevtALL/sumwtim) << " sec^-1\"," << std::endl
             << "\"NumMatrixElements(notNan)\": " << nevtALL - nnan << "," << std::endl
             << std::scientific
             << "\"MeanMatrixElemValue\": "
             << "\"" << std::to_string(meanelem) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"StdErrMatrixElemValue\": "
             << "\"" << std::to_string(stdelem/sqrt(nevtALL)) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"StdDevMatrixElemValue\": "
             << "\"" << std::to_string(stdelem)
             << " GeV^" << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"MinMatrixElemValue\": "
             << "\"" << std::to_string(minelem) << " GeV^"
             << std::to_string(meGeVexponent) << "\"," << std::endl
             << "\"MaxMatrixElemValue\": "
             << "\"" << std::to_string(maxelem) << " GeV^"
             << std::to_string(meGeVexponent) <<  "\"," << std::endl;

    timermap.dump(jsonFile, true); // NB For the active json timer this dumps a partial total

    jsonFile << "}" << std::endl << "]";
    jsonFile.close();
  }

  // *** STOP THE NEW TIMERS ***
  timermap.stop();
  if (perf)
  {
    std::cout << "***********************************************************************" << std::endl;
    timermap.dump();
    std::cout << "***********************************************************************" << std::endl;
  }

  //std::cout << "ALL OK" << std::endl;
}
