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

#ifndef ALPAKA
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#if defined(ALPAKA) || defined(__CUDACC__)
#include "rambo.cc"
#else
#include "rambo.h"
#endif

#ifdef MGONGPU_COMMONRAND_ONHOST
#include "CommonRandomNumbers.h"
#endif
#include "CPPProcess.h"
#include "Memory.h"
#include "timermap.h"

#include "epoch_process_id.h"
#define STRINGIFY(s) #s
#define XSTRINGIFY(s) STRINGIFY(s)

#define SEP79 79

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return (int)strlen(s) == t - s;
}

// Disabling fast math is essential here, otherwise results are undefined
// See https://stackoverflow.com/a/40702790 about __attribute__ on gcc
// See https://stackoverflow.com/a/32292725 about __attribute__ on clang
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
bool fp_is_abnormal( const fptype& fp )
{
  if ( std::isnan( fp ) ) return true;
  if ( fp != fp ) return true;
  return false;
}

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
bool fp_is_zero( const fptype& fp )
{
  if ( fp == 0 ) return true;
  return false;
}

// See https://en.cppreference.com/w/cpp/numeric/math/FP_categories
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
const char* fp_show_class( const fptype& fp )
{
  switch( std::fpclassify( fp ) ) {
  case FP_INFINITE:  return "Inf";
  case FP_NAN:       return "NaN";
  case FP_NORMAL:    return "normal";
  case FP_SUBNORMAL: return "subnormal";
  case FP_ZERO:      return "zero";
  default:           return "unknown";
  }
}

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
void debug_me_is_abnormal( const fptype& me, int ievtALL )
{
  std::cout << "DEBUG[" << ievtALL << "]"
            << " ME=" << me
            << " fpisabnormal=" << fp_is_abnormal( me )
            << " fpclass=" << fp_show_class( me )
            << " (me==me)=" << ( me == me )
            << " (me==me+1)=" << ( me == me+1 )
            << " isnan=" << std::isnan( me )
            << " isfinite=" << std::isfinite( me )
            << " isnormal=" << std::isnormal( me )
            << " is0=" << ( me == 0 )
            << " is1=" << ( me == 1 )
            << " abs(ME)=" << std::abs( me )
            << " isnan=" << std::isnan( std::abs( me ) )
            << std::endl;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j] [--json_file <JSON_FILE>] [--param_card <PARAM_CARD_FILE>]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations [jsondate] [jsonrun]" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only enables NaN/abnormal warnings and OMP debugging" << std::endl;
  std::cout << "The '--json_file' argument sets the path to the file where the JSON output is saved if '--json | -j' is specified. (default: ./perf/data/<jsondate>-perf-test-run<jsonrun>.json)" << std::endl; 
  std::cout << "The '--param_card' argument sets the path to the param_card file (default: ../../Cards/param_card.dat)" << std::endl;
#if !defined(ALPAKA) && !defined(__CUDACC__)
#ifdef _OPENMP
  std::cout << std::endl << "Use the OMP_NUM_THREADS environment variable to control OMP multi-threading" << std::endl;
  std::cout << "(OMP multithreading will be disabled if OMP_NUM_THREADS is not set)" << std::endl;
#endif
#endif
  std::cout << "The '--help|-h' flag prints this message" << std::endl;
  return ret;
}

int main(int argc, char **argv)
{
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
  std::string param_card = "../../Cards/param_card.dat";
  bool json_file_bool = false;
  std::string json_file = "";

  for (int argn = 1; argn < argc; ++argn) {
    std::string arg = argv[argn];
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
    else if ( ( arg == "--help" ) || ( arg == "-h" ) )
    {
      return usage(argv[0]);
    }
    else if ( arg == "--param_card" )
    {
      param_card = argv[argn + 1];
      argn++;
    }
    else if ( arg == "--json_file" )
    {
      json_file = argv[argn + 1];
      json_file_bool = true;
      argn++;
    }
    else if ( arg == "--device_id" )
    {
      //Do nothing, skip input
      argn++;
    }
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

  const int neppR = mgOnGpu::neppR; // AOSOA layout: constant at compile-time
  if ( gputhreads%neppR != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of neppR=" << neppR << std::endl;
    return usage(argv[0]);
  }

  const int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
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

  if ( gputhreads < 1 )
  {
    std::cout << "ERROR! #threads/block should be >= 0" << std::endl;
    return usage(argv[0]);
  }

  if ( gpublocks < 1 )
  {
    std::cout << "ERROR! #blocks should be >= 0" << std::endl;
    return usage(argv[0]);
  }

#if !defined(ALPAKA) && !defined(__CUDACC__)
#ifdef _OPENMP
  // Set OMP_NUM_THREADS equal to 1 if it is not yet set
  char* ompnthr = getenv( "OMP_NUM_THREADS" );
  if ( debug )
  {
    std::cout << "DEBUG: omp_get_num_threads() = " << omp_get_num_threads() << std::endl; // always == 1 here!
    std::cout << "DEBUG: omp_get_max_threads() = " << omp_get_max_threads() << std::endl;
    std::cout << "DEBUG: ${OMP_NUM_THREADS}    = '" << ( ompnthr == 0 ? "[not set]" : ompnthr ) << "'" << std::endl;
  }
  if ( ompnthr == NULL || std::string(ompnthr).find_first_not_of("0123456789") != std::string::npos || atol( ompnthr ) == 0 )
  {
    if ( ompnthr != NULL ) std::cout << "WARNING! OMP_NUM_THREADS is invalid: will use only 1 thread" << std::endl;
    else if ( debug ) std::cout << "DEBUG: OMP_NUM_THREADS is not set: will use only 1 thread" << std::endl;
    omp_set_num_threads( 1 ); // https://stackoverflow.com/a/22816325
    if ( debug ) std::cout << "DEBUG: omp_get_num_threads() = " << omp_get_num_threads() << std::endl; // always == 1 here!
    if ( debug ) std::cout << "DEBUG: omp_get_max_threads() = " << omp_get_max_threads() << std::endl;
  }

  // Fail gently and avoid "Illegal instruction (core dumped)" if the host does not support the requested AVX
  // [NB: this prevents a crash on pmpe04 but not on some github CI nodes]
  auto supportsAvx = [](){
#if defined __AVX512VL__
    bool ok = __builtin_cpu_supports( "avx512vl" );
    const std::string tag = "skylake-avx512 (AVX512VL)";
#elif defined __AVX2__
    bool ok = __builtin_cpu_supports( "avx2" );
    const std::string tag = "haswell (AVX2)";
#elif defined __SSE4_2__
#ifdef __PPC__
    // See https://gcc.gnu.org/onlinedocs/gcc/Basic-PowerPC-Built-in-Functions-Available-on-all-Configurations.html
    bool ok = __builtin_cpu_supports( "vsx" );
    const std::string tag = "powerpc vsx (128bit as in SSE4.2)";
#else
    bool ok = __builtin_cpu_supports( "sse4.2" );
    const std::string tag = "nehalem (SSE4.2)";
#endif
#else
    bool ok = true;
    const std::string tag = "none";
#endif
    if ( tag == "none" )
      std::cout << "INFO: The application does not require the host to support any AVX feature" << std::endl;
    else if ( ok )
      std::cout << "INFO: The application is built for " << tag << " and the host supports it" << std::endl;
    else
      std::cout << "ERROR! The application is built for " << tag << " but the host does not support it" << std::endl;
    return ok;
  };
  if ( ! supportsAvx() ) return 1;
#endif
#endif

  const int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  const int nevt = ndim; // number of events in one iteration == number of GPU threads
  const int nevtALL = niter*nevt; // total number of ALL events in all iterations

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // *** START THE NEW TIMERS ***
  mgOnGpu::TimerMap timermap;

  // === STEP 0 - INITIALISE

#if defined(ALPAKA) || defined(__CUDACC__)
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
      //if ( _print ) std::cout << "Calling cudaDeviceReset()." << std::endl;
      checkCuda( cudaDeviceReset() ); // this is needed by cuda-memcheck --leak-check full
    }
    bool _print{false};
  } cudaTearDown(debug);
#endif

  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
#if defined(ALPAKA) || defined(__CUDACC__)
  gProc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#else
  Proc::CPPProcess process( niter, gpublocks, gputhreads, verbose );
#endif

  // Read param_card and set parameters
  process.initProc(param_card);
  const fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
  //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
  const int meGeVexponent = -(2 * mgOnGpu::npar - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  using mgOnGpu::np4;
  using mgOnGpu::nparf;
  using mgOnGpu::npar;
  using mgOnGpu::ncomb; // Number of helicity combinations
  const int nRnarray = np4*nparf*nevt; // (NB: AOSOA layout with nevt=npagR*neppR events per iteration)
  const int nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for AOSOA layouts)
  const int nWeights = nevt;
  const int nMEs     = nevt; // FIXME: assume process.nprocesses == 1 (eventually: nMEs = nevt * nprocesses?)

#if defined MGONGPU_CURAND_ONHOST or defined MGONGPU_COMMONRAND_ONHOST or not (defined __CUDACC__ or defined ALPAKA)
  auto hstRnarray   = hstMakeUnique<fptype   >( nRnarray ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
#endif
  auto hstMomenta   = hstMakeUnique<fptype_sv>( nMomenta ); // AOSOA[npagM][npar][np4][neppM] (NB: nevt=npagM*neppM)
  auto hstIsGoodHel = hstMakeUnique<bool     >( ncomb );
  auto hstWeights   = hstMakeUnique<fptype   >( nWeights );
  auto hstMEs       = hstMakeUnique<fptype_sv>( nMEs ); // AOSOA[npagM][neppM] (NB: nevt=npagM*neppM)

#if defined(ALPAKA) || defined(__CUDACC__)
  auto devRnarray   = devMakeUnique<fptype   >( nRnarray ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  auto devMomenta   = devMakeUnique<fptype   >( nMomenta ); // AOSOA[npagM][npar][np4][neppM] (NB: nevt=npagM*neppM)
  auto devIsGoodHel = devMakeUnique<bool     >( ncomb );
  auto devWeights   = devMakeUnique<fptype   >( nWeights );
  auto devMEs       = devMakeUnique<fptype   >( nMEs ); // AOSOA[npagM][neppM] (NB: nevt=npagM*neppM)

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
  std::unique_ptr<double[]> wv3atimes( new double[niter] );
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
#if defined(ALPAKA) || defined(__CUDACC__)
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
#if defined(ALPAKA) || defined(__CUDACC__)
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
#elif defined(ALPAKA) || defined(__CUDACC__)
#ifdef MGONGPU_CURAND_ONDEVICE
    grambo2toNm0::generateRnarray( rnGen, devRnarray.get(), nevt );
#elif defined MGONGPU_CURAND_ONHOST
    grambo2toNm0::generateRnarray( rnGen, hstRnarray.get(), nevt );
#endif
#else
    rambo2toNm0::generateRnarray( rnGen, hstRnarray.get(), nevt );
#endif
    //std::cout << "Got random numbers" << std::endl;

#if defined(ALPAKA) || defined(__CUDACC__)
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
#ifdef ALPAKA
    CUPLA_KERNEL(grambo2toNm0::getMomentaInitial)(gpublocks, gputhreads)( energy, devMomenta.get() );
#elif defined __CUDACC__
    grambo2toNm0::getMomentaInitial<<<gpublocks, gputhreads>>>( energy, devMomenta.get() );
#else
    rambo2toNm0::getMomentaInitial( energy, hstMomenta.get(), nevt );
#endif
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );
#ifdef ALPAKA
    CUPLA_KERNEL(grambo2toNm0::getMomentaFinal)(gpublocks, gputhreads)( energy, devRnarray.get(), devMomenta.get(), devWeights.get() );
#elif defined __CUDACC__
    grambo2toNm0::getMomentaFinal<<<gpublocks, gputhreads>>>( energy, devRnarray.get(), devMomenta.get(), devWeights.get() );
#else
    rambo2toNm0::getMomentaFinal( energy, hstRnarray.get(), hstMomenta.get(), hstWeights.get(), nevt );
#endif
    //std::cout << "Got final momenta" << std::endl;

#if defined(ALPAKA) || defined(__CUDACC__)
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
    if ( iiter == 0 )
    {
      const std::string ghelKey = "0d SGoodHel";
      timermap.start( ghelKey );
#if defined(ALPAKA) || defined(__CUDACC__)
      // ... 0d1. Compute good helicity mask on the device
#ifdef ALPAKA
      CUPLA_KERNEL(gProc::sigmaKin_getGoodHel)(gpublocks, gputhreads)(devMomenta.get(), devMEs.get(), devIsGoodHel.get());
#else
      gProc::sigmaKin_getGoodHel<<<gpublocks, gputhreads>>>(devMomenta.get(), devMEs.get(), devIsGoodHel.get());
#endif
      checkCuda( cudaPeekAtLastError() );
      // ... 0d2. Copy back good helicity mask to the host
      checkCuda( cudaMemcpy( hstIsGoodHel.get(), devIsGoodHel.get(), nbytesIsGoodHel, cudaMemcpyDeviceToHost ) );
      // ... 0d3. Copy back good helicity list to constant memory on the device
      gProc::sigmaKin_setGoodHel(hstIsGoodHel.get());
#else
      // ... 0d1. Compute good helicity mask on the host
      Proc::sigmaKin_getGoodHel(hstMomenta.get(), hstMEs.get(), hstIsGoodHel.get(), nevt);
      // ... 0d2. Copy back good helicity list to static memory on the host
      Proc::sigmaKin_setGoodHel(hstIsGoodHel.get());
#endif
    }

    // *** START THE OLD-STYLE TIMERS FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    double wavetime = 0; // calc plus copy
    double wv3atime = 0; // calc only

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
#if defined(ALPAKA) || defined(__CUDACC__)
#ifndef MGONGPU_NSIGHT_DEBUG
#ifdef ALPAKA
    CUPLA_KERNEL(gProc::sigmaKin)(gpublocks, gputhreads)(devMomenta.get(), devMEs.get());
#else
    gProc::sigmaKin<<<gpublocks, gputhreads>>>(devMomenta.get(), devMEs.get());
#endif
#else
#ifdef ALPAKA
    CUPLA_KERNEL(gProc::sigmaKin)(gpublocks, gputhreads, ntpbMAX*sizeof(float))(devMomenta.get(), devMEs.get());
#else
    gProc::sigmaKin<<<gpublocks, gputhreads, ntpbMAX*sizeof(float)>>>(devMomenta.get(), devMEs.get());
#endif
#endif
    checkCuda( cudaPeekAtLastError() );
    checkCuda( cudaDeviceSynchronize() );
#else
    Proc::sigmaKin(hstMomenta.get(), hstMEs.get(), nevt);
#endif

    // *** STOP THE NEW OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wv3atime += timermap.stop(); // calc only
    wavetime += wv3atime; // calc plus copy

#if defined(ALPAKA) || defined(__CUDACC__)
    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    timermap.start( cmesKey );
    checkCuda( cudaMemcpy( hstMEs.get(), devMEs.get(), nbytesMEs, cudaMemcpyDeviceToHost ) );
    // *** STOP THE OLD OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wavetime += timermap.stop(); // calc plus copy
#endif

    // === STEP 4 FINALISE LOOP
    // --- 4a Dump within the loop
    const std::string loopKey = "4a DumpLoop";
    timermap.start(loopKey);
    genrtimes[iiter] = genrtime;
    rambtimes[iiter] = rambtime;
    wavetimes[iiter] = wavetime;
    wv3atimes[iiter] = wv3atime;

    if (verbose)
    {
      std::cout << std::string(SEP79, '*') << std::endl
                << "Iteration #" << iiter+1 << " of " << niter << std::endl;
      if (perf) std::cout << "Wave function time: " << wavetime << std::endl;
    }

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
#ifndef MGONGPU_CPPSIMD
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*np4*neppM + 0*neppM + ieppM] // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*np4*neppM + 1*neppM + ieppM] // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*np4*neppM + 2*neppM + ieppM] // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*np4*neppM + 3*neppM + ieppM] // AOSOA[ipagM][ipar][3][ieppM]
#else
                    << std::setw(14) << hstMomenta[ipagM*npar*np4 + ipar*np4 + 0][ieppM] // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4 + ipar*np4 + 1][ieppM] // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4 + ipar*np4 + 2][ieppM] // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4 + ipar*np4 + 3][ieppM] // AOSOA[ipagM][ipar][3][ieppM]
#endif
                    << std::endl
                    << std::defaultfloat; // default format: affects all floats
        }
        std::cout << std::string(SEP79, '-') << std::endl;
        // Display matrix elements
        std::cout << " Matrix element = "
#ifndef MGONGPU_CPPSIMD
                  << hstMEs[ievt]
#else
                  << hstMEs[ievt/neppM][ievt%neppM]
#endif
                  << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(SEP79, '-') << std::endl;
      }
      // Fill the arrays with ALL MEs and weights
#ifndef MGONGPU_CPPSIMD
      matrixelementALL[iiter*nevt + ievt] = hstMEs[ievt]; // FIXME: assume process.nprocesses == 1
#else
      matrixelementALL[iiter*nevt + ievt] = hstMEs[ievt/neppM][ievt%neppM]; // FIXME: assume process.nprocesses == 1
#endif
      weightALL[iiter*nevt + ievt] = hstWeights[ievt];
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

  double sumw3atim = 0;
  double sqsw3atim = 0;
  double minw3atim = wv3atimes[0];
  double maxw3atim = wv3atimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumw3atim += wv3atimes[iiter];
    sqsw3atim += wv3atimes[iiter]*wv3atimes[iiter];
    minw3atim = std::min( minw3atim, wv3atimes[iiter] );
    maxw3atim = std::max( maxw3atim, wv3atimes[iiter] );
  }
  double meanw3atim = sumw3atim / niter;
  //double stdw3atim = std::sqrt( sqsw3atim / niter - meanw3atim * meanw3atim );

  int nabn = 0;
  int nzero = 0;
  double minelem = matrixelementALL[0];
  double maxelem = matrixelementALL[0];
  double minweig = weightALL[0];
  double maxweig = weightALL[0];
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // The following events are abnormal in a run with "-p 2048 256 12 -d"
    // - check.exe/commonrand: ME[310744,451171,3007871,3163868,4471038,5473927] with fast math
    // - check.exe/curand: ME[578162,1725762,2163579,5407629,5435532,6014690] with fast math
    // - gcheck.exe/curand: ME[596016,1446938] with fast math
    // Debug NaN/abnormal issues
    //if ( ievtALL == 310744 ) // this ME is abnormal both with and without fast math
    //  debug_me_is_abnormal( matrixelementALL[ievtALL], ievtALL );
    //if ( ievtALL == 5473927 ) // this ME is abnormal only with fast math
    //  debug_me_is_abnormal( matrixelementALL[ievtALL], ievtALL );
    // Compute min/max
    if ( fp_is_zero( matrixelementALL[ievtALL] ) ) nzero++;
    if ( fp_is_abnormal( matrixelementALL[ievtALL] ) )
    {
      if ( debug ) // only printed out with "-p -d" (matrixelementALL is not filled without -p)
        std::cout << "WARNING! ME[" << ievtALL << "] is NaN/abnormal" << std::endl;
      nabn++;
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
    if ( fp_is_abnormal( matrixelementALL[ievtALL] ) ) continue;
    sumelemdiff += ( matrixelementALL[ievtALL] - minelem );
    sumweigdiff += ( weightALL[ievtALL] - minweig );
  }
  double meanelem = minelem + sumelemdiff / ( nevtALL - nabn );
  double meanweig = minweig + sumweigdiff / ( nevtALL - nabn );
  double sqselemdiff = 0;
  double sqsweigdiff = 0;
  for ( int ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
    // Compute stddev from the squared sum of diff to mean
    if ( fp_is_abnormal( matrixelementALL[ievtALL] ) ) continue;
    sqselemdiff += std::pow( matrixelementALL[ievtALL] - meanelem, 2 );
    sqsweigdiff += std::pow( weightALL[ievtALL] - meanweig, 2 );
  }
  double stdelem = std::sqrt( sqselemdiff / ( nevtALL - nabn ) );
  double stdweig = std::sqrt( sqsweigdiff / ( nevtALL - nabn ) );

  // === STEP 9 FINALISE
  // --- 9a. Destroy curand generator
  const std::string dgenKey = "9a GenDestr";
  timermap.start( dgenKey );
#ifndef MGONGPU_COMMONRAND_ONHOST
#if defined(ALPAKA) || defined(__CUDACC__)
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
#if !defined(ALPAKA) && !defined(__CUDACC__)
    // Get the output of "nproc --all" (https://stackoverflow.com/a/478960)
    std::string nprocall;
    std::unique_ptr<FILE, decltype(&pclose)> nprocpipe( popen( "nproc --all", "r" ), pclose );
    if ( !nprocpipe ) throw std::runtime_error( "`nproc --all` failed?" );
    std::array<char, 128> nprocbuf;
    while ( fgets( nprocbuf.data(), nprocbuf.size(), nprocpipe.get() ) != nullptr ) nprocall += nprocbuf.data();
#endif
#ifdef MGONGPU_CPPSIMD
#ifdef MGONGPU_HAS_CXTYPE_REF
    const std::string cxtref = " [cxtype_ref=YES]";
#else
    const std::string cxtref = " [cxtype_ref=NO]";
#endif
#endif
    // Dump all configuration parameters and all results
    std::cout << std::string(SEP79, '*') << std::endl
#ifdef ALPAKA
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_ALPAKA"
#elif defined __CUDACC__
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_CUDA"
#else
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_CPP"
#endif
              << " [" << process.getCompiler() << "]"
#ifdef MGONGPU_INLINE_HELAMPS
              << " [inlineHel=1]" << std::endl
#else
              << " [inlineHel=0]" << std::endl
#endif
              << "NumBlocksPerGrid            = " << gpublocks << std::endl
              << "NumThreadsPerBlock          = " << gputhreads << std::endl
              << "NumIterations               = " << niter << std::endl
              << std::string(SEP79, '-') << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
              << "FP precision                = DOUBLE (NaN/abnormal=" << nabn << ", zero=" << nzero << ")" << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
              << "FP precision                = FLOAT (NaN/abnormal=" << nabn << ", zero=" << nzero << ")" << std::endl
#endif
#if defined(ALPAKA) || defined(__CUDACC__)
#if defined MGONGPU_CXTYPE_CUCOMPLEX
              << "Complex type                = CUCOMPLEX" << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
              << "Complex type                = THRUST::COMPLEX" << std::endl
#elif defined MGONGPU_CXTYPE_ALSIMPLE
              << "Complex type                = ALSIMPLE::COMPLEX" << std::endl
#endif
#else
              << "Complex type                = STD::COMPLEX" << std::endl
#endif
              << "RanNumb memory layout       = AOSOA[" << neppR << "]"
              << ( neppR == 1 ? " == AOS" : "" )
              << " [HARDCODED FOR REPRODUCIBILITY]" << std::endl
              << "Momenta memory layout       = AOSOA[" << neppM << "]"
              << ( neppM == 1 ? " == AOS" : "" ) << std::endl
#if defined(ALPAKA) || defined(__CUDACC__)
      //<< "Wavefunction GPU memory     = LOCAL" << std::endl
#else
#if !defined MGONGPU_CPPSIMD
              << "Internal loops fptype_sv    = SCALAR ('none': ~vector[" << neppV
              << "], no SIMD)" << std::endl
#elif defined __AVX512VL__
#ifdef MGONGPU_PVW512
              << "Internal loops fptype_sv    = VECTOR[" << neppV
              << "] ('512z': AVX512, 512bit)" << cxtref << std::endl
#else
              << "Internal loops fptype_sv    = VECTOR[" << neppV
              << "] ('512y': AVX512, 256bit)" << cxtref << std::endl
#endif
#elif defined __AVX2__
              << "Internal loops fptype_sv    = VECTOR[" << neppV
              << "] ('avx2': AVX2, 256bit)" << cxtref << std::endl
#elif defined __SSE4_2__
              << "Internal loops fptype_sv    = VECTOR[" << neppV
#ifdef __PPC__
              << "] ('sse4': PPC VSX, 128bit)" << cxtref << std::endl
#else
              << "] ('sse4': SSE4.2, 128bit)" << cxtref << std::endl
#endif
#else
#error Internal error: unknown SIMD build configuration
#endif
#endif
#ifdef ALPAKA
#if defined MGONGPU_COMMONRAND_ONHOST
              << "Random number generation    = COMMON RANDOM HOST (ALPAKA code)" << std::endl
#elif defined MGONGPU_CURAND_ONDEVICE
              << "Random number generation    = CURAND DEVICE (ALPAKA code)" << std::endl
#elif defined MGONGPU_CURAND_ONHOST
              << "Random number generation    = CURAND HOST (ALPAKA code)" << std::endl
#endif
#elif defined(__CUDACC__)
#if defined MGONGPU_COMMONRAND_ONHOST
              << "Random number generation    = COMMON RANDOM HOST (CUDA code)" << std::endl
#elif defined MGONGPU_CURAND_ONDEVICE
              << "Random number generation    = CURAND DEVICE (CUDA code)" << std::endl
#elif defined MGONGPU_CURAND_ONHOST
              << "Random number generation    = CURAND HOST (CUDA code)" << std::endl
#endif
#else
#if defined MGONGPU_COMMONRAND_ONHOST
              << "Random number generation    = COMMON RANDOM (C++ code)" << std::endl
#else
              << "Random number generation    = CURAND (C++ code)" << std::endl
#endif
#ifdef _OPENMP
              << "OMP threads / `nproc --all` = " << omp_get_max_threads() << " / " << nprocall // includes a newline
#endif
#endif
      //<< "MatrixElements compiler     = " << process.getCompiler() << std::endl
              << std::string(SEP79, '-') << std::endl
              << "NumberOfEntries             = " << niter << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "TotalTime[Rnd+Rmb+ME] (123) = ( " << sumgtim+sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo+ME]    (23) = ( " << sumrtim+sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[RndNumGen]    (1) = ( " << sumgtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[Rambo]        (2) = ( " << sumrtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[MatrixElems]  (3) = ( " << sumwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "MeanTimeInMatrixElems       = ( " << meanwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "[Min,Max]TimeInMatrixElems  = [ " << minwtim
              << " ,  " << maxwtim << " ]  sec" << std::endl
      //<< "StdDevTimeInMatrixElems     = ( " << stdwtim << std::string(16, ' ') << " )  sec" << std::endl
              << "TotalTime[MECalcOnly]  (3a) = ( " << sumw3atim << std::string(16, ' ') << " )  sec" << std::endl
              << "MeanTimeInMECalcOnly        = ( " << meanw3atim << std::string(16, ' ') << " )  sec" << std::endl
              << "[Min,Max]TimeInMECalcOnly   = [ " << minw3atim
              << " ,  " << maxw3atim << " ]  sec" << std::endl
      //<< "StdDevTimeInMECalcOnly      = ( " << stdw3atim << std::string(16, ' ') << " )  sec" << std::endl
              << std::string(SEP79, '-') << std::endl
      //<< "ProcessID:                  = " << getpid() << std::endl
      //<< "NProcesses                  = " << process.nprocesses << std::endl
              << "TotalEventsComputed         = " << nevtALL << std::endl
              << "EvtsPerSec[Rnd+Rmb+ME](123) = ( " << nevtALL/(sumgtim+sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[Rmb+ME]     (23) = ( " << nevtALL/(sumrtim+sumwtim)
              << std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[RndNumbGen]   (1) = ( " << nevtALL/sumgtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
      //<< "EvtsPerSec[Rambo]        (2) = ( " << nevtALL/sumrtim
      //<< std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[MatrixElems] (3) = ( " << nevtALL/sumwtim
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << "EvtsPerSec[MECalcOnly] (3a) = ( " << nevtALL/sumw3atim
              << std::string(16, ' ') << " )  sec^-1" << std::endl
              << std::defaultfloat; // default format: affects all floats
    std::cout << std::string(SEP79, '*') << std::endl
              << "NumMatrixElems(notAbnormal) = " << nevtALL - nabn << std::endl
              << std::scientific // fixed format: affects all floats (default precision: 6)
              << "MeanMatrixElemValue         = ( " << meanelem
              << " +- " << stdelem/sqrt(nevtALL - nabn) << " )  GeV^" << meGeVexponent << std::endl // standard error
              << "[Min,Max]MatrixElemValue    = [ " << minelem
              << " ,  " << maxelem << " ]  GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue       = ( " << stdelem << std::string(16, ' ') << " )  GeV^" << meGeVexponent << std::endl
              << "MeanWeight                  = ( " << meanweig
              << " +- " << stdweig/sqrt(nevtALL - nabn) << " )" << std::endl // standard error
              << "[Min,Max]Weight             = [ " << minweig
              << " ,  " << maxweig << " ]" << std::endl
              << "StdDevWeight                = ( " << stdweig << std::string(16, ' ') << " )" << std::endl
              << std::defaultfloat; // default format: affects all floats
  }

  // --- 9c Dump to json
  const std::string jsonKey = "9c DumpJson";
  timermap.start(jsonKey);

  if(json)
  {
    std::string jsonFileName = std::to_string(jsondate) + "-perf-test-run" + std::to_string(jsonrun) + ".json";
    jsonFileName = "./perf/data/" + jsonFileName;
    if (json_file_bool) {
        jsonFileName = json_file;
    }

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
      if ( system(command) != 0 )
        std::cout << "WARNING! Command '" << temp << "' failed" << std::endl;
      jsonFile << ", " << std::endl;
    }

    jsonFile << "{" << std::endl
             << "\"NumIterations\": " << niter << ", " << std::endl
             << "\"NumThreadsPerBlock\": " << gputhreads << ", " << std::endl
             << "\"NumBlocksPerGrid\": " << gpublocks << ", " << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
             << "\"FP precision\": " << "\"DOUBLE (NaN/abnormal=" << nabn << ")\"," << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
             << "\"FP precision\": " << "\"FLOAT (NaN/abnormal=" << nabn << ")\"," << std::endl
#endif
             << "\"Complex type\": "
#if defined(ALPAKA) || defined(__CUDACC__)
#if defined MGONGPU_CXTYPE_CUCOMPLEX
             << "\"CUCOMPLEX\"," << std::endl
#elif defined MGONGPU_CXTYPE_THRUST
             << "\"THRUST::COMPLEX\"," << std::endl
#elif defined MGONGPU_CXTYPE_ALSIMPLE
             << "\"ALSIMPLE::COMPLEX\"," << std::endl
#endif
#else
             << "\"STD::COMPLEX\"," << std::endl
#endif
             << "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\""
             << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Momenta memory layout\": " << "\"AOSOA[" << neppM << "]\""
             << ( neppM == 1 ? " == AOS" : "" ) << ", " << std::endl
#if defined(ALPAKA) || defined(__CUDACC__)
      //<< "\"Wavefunction GPU memory\": " << "\"LOCAL\"," << std::endl
#endif
             << "\"Curand generation\": "
#ifdef ALPAKA
#if defined MGONGPU_COMMONRAND_ONHOST
             << "\"COMMON RANDOM HOST (ALPAKA code)\"," << std::endl;
#elif defined MGONGPU_CURAND_ONDEVICE
    << "\"CURAND DEVICE (ALPAKA code)\"," << std::endl;
#elif defined MGONGPU_CURAND_ONHOST
    << "\"CURAND HOST (ALPAKA code)\"," << std::endl;
#endif
#elif defined(__CUDACC__)
#if defined MGONGPU_COMMONRAND_ONHOST
             << "\"COMMON RANDOM HOST (CUDA code)\"," << std::endl;
#elif defined MGONGPU_CURAND_ONDEVICE
    << "\"CURAND DEVICE (CUDA code)\"," << std::endl;
#elif defined MGONGPU_CURAND_ONHOST
    << "\"CURAND HOST (CUDA code)\"," << std::endl;
#endif
#else
#if defined MGONGPU_COMMONRAND_ONHOST
    << "\"COMMON RANDOM (C++ code)\"," << std::endl;
#else
    << "\"CURAND (C++ code)\"," << std::endl;
#endif
#endif
    jsonFile << "\"NumberOfEntries\": " << niter << "," << std::endl
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
             << std::to_string(nevtALL/(sumgtim+sumrtim+sumwtim)) << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[Rmb+ME] (23)\": \""
             << std::to_string(nevtALL/(sumrtim+sumwtim)) << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[MatrixElems] (3)\": \""
             << std::to_string(nevtALL/sumwtim) << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[MECalcOnly] (3)\": \""
             << std::to_string(nevtALL/sumw3atim) << " sec^-1\"," << std::endl
             << "\"NumMatrixElems(notAbnormal)\": " << nevtALL - nabn << "," << std::endl
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
    std::cout << std::string(SEP79, '*') << std::endl;
    timermap.dump();
    std::cout << std::string(SEP79, '*') << std::endl;
  }

  //std::cout << "ALL OK" << std::endl;
  return 0;
}
