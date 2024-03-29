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

#include "Kokkos_Core.hpp"
#include "mgOnGpuTypes.h"
#include "random_generator.h"
#include "rambo.h"

#include "CPPProcess.h"
#include "timermap.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "epoch_process_id.h"
#define STRINGIFY(s) #s
#define XSTRINGIFY(s) STRINGIFY(s)

#define SEP79 79

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "Bridge.h"
#endif

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args ) {
    //See https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
    int size_s = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    auto buf = std::make_unique<char[]>( size );
    snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return (int)strlen(s) == t - s;
}

// Disabling fast math is essential here, otherwise results are undefined
// See https://stackoverflow.com/a/40702790 about __attribute__ on gcc
// See https://stackoverflow.com/a/32292725 about __attribute__ on clang
__attribute__((optimize("-fno-fast-math")))
bool fp_is_abnormal( const fptype& fp )
{
  if ( std::isnan( fp ) ) return true;
  if ( fp != fp ) return true;
  return false;
}

__attribute__((optimize("-fno-fast-math")))
bool fp_is_zero( const fptype& fp )
{
  if ( fp == 0 ) return true;
  return false;
}

// See https://en.cppreference.com/w/cpp/numeric/math/FP_categories
__attribute__((optimize("-fno-fast-math")))
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

__attribute__((optimize("-fno-fast-math")))
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
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j] [--json_file <JSON_FILE>] [--param_card <PARAM_CARD_FILE>] [--bridge] [--device_info] --device_id <DEVICE_ID>"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations [jsondate] [jsonrun]" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only enables NaN/abnormal warnings and OMP debugging" << std::endl;
  std::cout << "The '--json_file' argument sets the path to the file where the JSON output is saved if '--json | -j' is specified. (default: ./perf/data/<jsondate>-perf-test-run<jsonrun>.json)" << std::endl; 
  std::cout << "The '--param_card' argument sets the path to the param_card file (default: ../../Cards/param_card.dat)" << std::endl;
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
  size_t niter = 0;
  size_t league_size = 1;
  size_t team_size = 32;
  size_t jsondate = 0;
  size_t jsonrun = 0;
  size_t numvec[5] = {0,0,0,0,0};
  size_t nnum = 0;
  std::string param_card = "../../Cards/param_card.dat";
  bool json_file_bool = false;
  std::string json_file = "";
  bool bridge = false;

  for (size_t argn = 1; argn < argc; ++argn) {
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
    else if ( arg == "--bridge" )
    {
      bridge = true;
    }
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
    league_size = numvec[0];
    team_size = numvec[1];
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
  // initialize Kokkos
  Kokkos::initialize(argc, argv); {


  using mgOnGpu::ntpbMAX;
  if ( team_size > ntpbMAX )
  {
    std::cout << "ERROR! #threads/block should be <= " << ntpbMAX << std::endl;
    return usage(argv[0]);
  }

  if ( team_size < 1 )
  {
    std::cout << "ERROR! #threads/block should be >= 0" << std::endl;
    return usage(argv[0]);
  }

  if ( league_size < 1 )
  {
    std::cout << "ERROR! #blocks should be >= 0" << std::endl;
    return usage(argv[0]);
  }

  const size_t ndim = league_size * team_size; // number of threads in one GPU grid
  const size_t nevt = ndim; // number of events in one iteration == number of GPU threads
  const size_t nevtALL = niter*nevt; // total number of ALL events in all iterations

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // *** START THE NEW TIMERS ***
  mgOnGpu::TimerMap timermap;

  // === STEP 0 - INITIALISE
  // --- 00. Initialise cuda (call cudaFree to ease cuda profile analysis)
  
  // --- 00. Initialise cuda (call cudaFree to ease cuda profile analysis)
  const std::string cdfrKey = "00 CudaFree";
  timermap.start( cdfrKey );
  
  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
  CPPProcess process( verbose );


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
  // const int nRnarray = np4*nparf*nevt; // (NB: AOSOA layout with nevt=npagR*neppR events per iteration)
  // const int nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for AOSOA layouts)
  // const int nWeights = nevt;
  // const int nMEs     = nevt; // FIXME: assume process.nprocesses == 1 (eventually: nMEs = nevt * nprocesses?)
  
  // Random Numbers
  Kokkos::View<fptype**> devRnarray(Kokkos::ViewAllocateWithoutInitializing("devRnarray"),nevt,np4*nparf);
#ifdef FIXED_RANDOM
    auto hstRnarray = Kokkos::create_mirror_view(devRnarray);
#endif

  // momenta
  Kokkos::View<fptype***> devMomenta(Kokkos::ViewAllocateWithoutInitializing("devMomenta"),nevt,npar,np4);
  auto hstMomenta = Kokkos::create_mirror_view(devMomenta);

  // matrix elements
  Kokkos::View<fptype*> devMEs(Kokkos::ViewAllocateWithoutInitializing("devMEs"),nevt);
  auto hstMEs = Kokkos::create_mirror_view(devMEs);

  constexpr fptype fixedG = 1.2177157847767195; // fixed G for aS=0.118 (hardcoded for now in check_sa.cc, fcheck_sa.f, 
  Kokkos::View<fptype*> devGs(Kokkos::ViewAllocateWithoutInitializing("devGs"),nevt);
  auto hstGs = Kokkos::create_mirror_view(devGs);
  for(size_t i=0; i < nevt; i++) hstGs[i] = fixedG;
  Kokkos::deep_copy(devGs,hstGs);

  // Set couplings and parameters
  Kokkos::View<cxtype*> dev_independent_couplings(Kokkos::ViewAllocateWithoutInitializing("dev_independent_couplings"), independentCouplings::nicoup );
  Kokkos::View<fptype*> dev_independent_parameters(Kokkos::ViewAllocateWithoutInitializing("dev_independent_parameters"), mgOnGpu::nparams );

  #if MGONGPU_NICOUP > 0
      //Create Kokkos view of independent couplings on host
      Kokkos::View<cxtype*, Kokkos::HostSpace> hst_independent_couplings(process.get_tIPC_ptr(), MGONGPU_NICOUP);
      
      //FIXME deep_copy doesn't work since dev_independent_couplings and hst_independent_couplings are not in the same execution space
      //Trying this "workaround"
      //Create mirror view of dev_independent_couplings on host 
      typename Kokkos::View<cxtype*>::HostMirror hst_mirror_independent_couplings = Kokkos::create_mirror_view(dev_independent_couplings);

      //Copy the data from host to device
      Kokkos::deep_copy(hst_mirror_independent_couplings, hst_independent_couplings); //Same execution space, so can copy. Should result in a no-op?
      Kokkos::deep_copy(dev_independent_couplings, hst_mirror_independent_couplings); //Copy from mirror
  #endif

  //Create Kokkos view of independent parameters on host
  Kokkos::View<fptype*, Kokkos::HostSpace> hst_independent_parameters(process.get_tIPD_ptr(), mgOnGpu::nparams);
  
  //FIXME deep_copy doesn't work since dev_independent_parameters and hst_independent_parameters are not in the same execution space
  //Trying this "workaround"
  //Create mirror view of dev_independent_parameters on host 
  typename Kokkos::View<fptype*>::HostMirror hst_mirror_independent_parameters = Kokkos::create_mirror_view(dev_independent_parameters);

  //Copy the data from host to device
  Kokkos::deep_copy(hst_mirror_independent_parameters, hst_independent_parameters); //Same execution space, so can copy. Should result in a no-op?
  Kokkos::deep_copy(dev_independent_parameters, hst_mirror_independent_parameters); //Copy from mirror

  // weights
  Kokkos::View<fptype*> devWeights(Kokkos::ViewAllocateWithoutInitializing("devWeights"),nevt);
  auto hstWeights = Kokkos::create_mirror_view(devWeights);

  // good helicity indices tracking (device-side only)
  Kokkos::View<size_t*> devNGoodHel("devNGoodHel",1); // TODO Fixed to 1 process
  auto hstNGoodHel = Kokkos::create_mirror_view(devNGoodHel);
  Kokkos::View<size_t*> devGoodHel("devGoodHel",ncomb);

  std::unique_ptr<double[]> genrtimes( new double[niter] );
  std::unique_ptr<double[]> rambtimes( new double[niter] );
  std::unique_ptr<double[]> wavetimes( new double[niter] );
  std::unique_ptr<double[]> wv3atimes( new double[niter] );
  std::unique_ptr<fptype[]> matrixelementALL( new fptype[nevtALL] ); // FIXME: assume process.nprocesses == 1
  std::unique_ptr<fptype[]> weightALL( new fptype[nevtALL] );

  double bridge_time = 0.0;
  #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      mg5amcKokkos::Bridge<fptype>* pbridge;
      fptype* _gs;
      fptype* _mes;
      fptype* _momenta;
      if (bridge) {
          //Initialize bridge object
          pbridge = new mg5amcKokkos::Bridge<fptype>( nevt, npar, np4 );
      
          _gs      = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _mes     = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _momenta = reinterpret_cast<fptype*>(std::malloc(nevt*npar*np4*sizeof(fptype)));
          for (size_t i=0; i < nevt*npar*np4; i++) {
              _momenta[i] = hstMomenta.data()[i];
          }
          for (size_t i=0; i < nevt; i++) {
              _gs[i] = hstGs.data()[i];
          }
      }
  #endif

  // --- 0c. Create curand or common generator
  const std::string cgenKey = "0c GenCreat";
  timermap.start( cgenKey );

  // init random number generator pool
  int random_seed = 1234;
  auto rand_pool = init_random_generator(random_seed);

  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************
  
  for (size_t iiter = 0; iiter < niter; ++iiter)
  {
    //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;

    // === STEP 1 OF 3

    // *** START THE OLD-STYLE TIMER FOR RANDOM GEN ***
    double genrtime = 0;

    // --- 1a. Seed curand generator (to get same results on host and device)
    // [NB This should not be necessary using the host API: "Generation functions
    // can be called multiple times on the same generator to generate successive
    // blocks of results. For pseudorandom generators, multiple calls to generation
    // functions will yield the same result as a single call with a large size."]
    const std::string sgenKey = "1a GenSeed ";
    timermap.start( sgenKey );
    // const unsigned long long seed = 20200805;

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );
#ifdef FIXED_RANDOM
    srand(123);
    
    for(size_t ii = 0; ii < nevt; ii++) {
        for(size_t jj = 0; jj < nparf; jj++) {
            for(size_t kk = 0; kk < np4; kk++) {
                fptype r = static_cast<double>(rand()) / (double)RAND_MAX;
                hstRnarray(ii, jj*np4 + kk) = r;
    }}}
    
    Kokkos::deep_copy(devRnarray,hstRnarray);
#else
    fill_random_numbers_2d(devRnarray,nevt,np4*nparf, rand_pool, league_size, team_size);
    Kokkos::fence();
#endif
    //std::cout << "Got random numbers" << std::endl;

    // *** STOP THE OLD-STYLE TIMER FOR RANDOM GEN ***
    genrtime += timermap.stop();

    // === STEP 2 OF 3
    // Fill in particle momenta for each of nevt events on the device

    // *** START THE OLD-STYLE TIMER FOR RAMBO ***
    double rambtime = 0;

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );
    get_initial_momenta(devMomenta,mgOnGpu::npar,energy,league_size,team_size);
    Kokkos::fence();
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );
    get_final_momenta(mgOnGpu::npari, mgOnGpu::npar, energy, devMomenta, devRnarray, devWeights, league_size, team_size);
    Kokkos::fence();
    //std::cout << "Got final momenta" << std::endl;

    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    rambtime += timermap.start( cwgtKey );
    Kokkos::deep_copy(hstWeights,devWeights);

    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    rambtime += timermap.start( cmomKey );
    Kokkos::deep_copy(hstMomenta,devMomenta);

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
      // ... 0d1. Compute good helicity mask on the device
      sigmaKin_setup(devMomenta, devGoodHel, devNGoodHel, devGs, dev_independent_couplings, dev_independent_parameters, league_size, team_size);
      Kokkos::fence();
    }

    // *** START THE OLD-STYLE TIMERS FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    double wavetime = 0; // calc plus copy
    double wv3atime = 0; // calc only
    
    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
    #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        sigmaKin(devMomenta, 0, devGoodHel, devNGoodHel, devGs, dev_independent_couplings, dev_independent_parameters, league_size, team_size, devMEs);
    #else
        sigmaKin(devMomenta, devGoodHel, devNGoodHel, devGs, dev_independent_couplings, dev_independent_parameters, league_size, team_size, devMEs);
    #endif
    Kokkos::fence();
    // *** STOP THE NEW OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wv3atime += timermap.stop(); // calc only
    wavetime += wv3atime; // calc plus copy

    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    timermap.start( cmesKey );
    Kokkos::deep_copy(hstMEs,devMEs);
    // *** STOP THE OLD OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wavetime += timermap.stop(); // calc plus copy

    #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        if (bridge) {
            const std::string bridgeKey = "3c bridge";
            timermap.start( bridgeKey );
            pbridge->gpu_sequence( _momenta, _gs, _mes, 0 );
            bridge_time += timermap.stop();
            
            //std::cout << "SA_MEs    Bridge_MEs    SA/Bridge" << std::endl;
            //for (size_t i=0; i < nevt; i++) {
            //    std::cout << hstMEs(i) << "    " << _mes[i] << std::endl;
            //}
        }
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

    for (size_t ievt = 0; ievt < nevt; ++ievt) // Loop over all events in this iteration
    {
      if (verbose)
      {
        // Display momenta
        std::cout << "Momenta:" << std::endl;
        for (size_t ipar = 0; ipar < npar; ipar++)
        {
          // NB: 'setw' affects only the next field (of any type)
          std::cout << std::scientific // fixed format: affects all floats (default precision: 6)
                    << std::setw(4) << ipar + 1
                    << std::setw(14) << hstMomenta(ievt,ipar,0) // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << hstMomenta(ievt,ipar,1) // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << hstMomenta(ievt,ipar,2) // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << hstMomenta(ievt,ipar,3) // AOSOA[ipagM][ipar][3][ieppM]
                    << std::endl
                    << std::defaultfloat; // default format: affects all floats
        }
        std::cout << std::string(SEP79, '-') << std::endl;
        // Display matrix elements
        std::cout << " Matrix element = "
                  << hstMEs[ievt]
                  << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(SEP79, '-') << std::endl;
      }
      // Fill the arrays with ALL MEs and weights
      matrixelementALL[iiter*nevt + ievt] = hstMEs[ievt]; // FIXME: assume process.nprocesses == 1
      weightALL[iiter*nevt + ievt] = hstWeights[ievt];
    }

    if (!(verbose || debug || perf))
    {
      std::cout << ".";
    }

  } // end loop iterations

  #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if (bridge) {
          // Destroy bridge objects
          std::free(_gs);
          std::free(_mes);
          std::free(_momenta);
          delete pbridge;
      }
  #endif

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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ievtALL++ )
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ievtALL++ )
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ievtALL++ )
  {
    // Compute stddev from the squared sum of diff to mean
    if ( fp_is_abnormal( matrixelementALL[ievtALL] ) ) continue;
    sqselemdiff += std::pow( matrixelementALL[ievtALL] - meanelem, 2 );
    sqsweigdiff += std::pow( weightALL[ievtALL] - meanweig, 2 );
  }
  double stdelem = std::sqrt( sqselemdiff / ( nevtALL - nabn ) );
  double stdweig = std::sqrt( sqsweigdiff / ( nevtALL - nabn ) );
  double bridge_MEs_per_sec = 0.0;
  if (bridge) {
      bridge_MEs_per_sec = nevtALL/bridge_time;
  }

  // === STEP 9 FINALISE

  Kokkos::deep_copy(hstNGoodHel,devNGoodHel);

  // --- 9a. Destroy curand generator
  const std::string dgenKey = "9a GenDestr";
  timermap.start( dgenKey );

  // --- 9b Dump to screen
  const std::string dumpKey = "9b DumpScrn";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
    // Dump all configuration parameters and all results
    std::cout << std::string(SEP79, '*') << std::endl
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_Kokkos"

#ifdef MGONGPU_INLINE_HELAMPS
              << " [inlineHel=1]" << std::endl
#else
              << " [inlineHel=0]" << std::endl
#endif
              << "NumBlocksPerGrid            = " << league_size << std::endl
              << "NumThreadsPerBlock          = " << team_size << std::endl
              << "NumIterations               = " << niter << std::endl
              << "IsGoodHel                   = " << hstNGoodHel(0) << " out of " << ncomb << std::endl
              << std::string(SEP79, '-') << std::endl
              << "FP precision                = " << FPTYPE_NAME << " (NaN/abnormal=" << nabn << ", zero=" << nzero << ")" << std::endl
              << "Complex type                = " << COMPLEX_TYPE_NAME << std::endl

              //<< "MatrixElements compiler     = " << process.getCompiler() << std::endl
              << std::string(SEP79, '-') << std::endl
              << "NumberOfEntries             = " << niter << std::endl
              << "TotalTime[Rnd+Rmb+ME] (123) = ( " << string_format("%1.15e", sumgtim+sumrtim+sumwtim) << " )  sec" << std::endl
              << "TotalTime[Rambo+ME]    (23) = ( " << string_format("%1.15e", sumrtim+sumwtim)         << " )  sec" << std::endl
              << "TotalTime[RndNumGen]    (1) = ( " << string_format("%1.15e", sumgtim)                 << " )  sec" << std::endl
              << "TotalTime[Rambo]        (2) = ( " << string_format("%1.15e", sumrtim)                 << " )  sec" << std::endl
              << "TotalTime[MatrixElems]  (3) = ( " << string_format("%1.15e", sumwtim)                 << " )  sec" << std::endl
              << "TotalTime[Bridge]           = ( " << string_format("%1.15e", bridge_time)             << " )  sec" << std::endl
              << "MeanTimeInMatrixElems       = ( " << string_format("%1.15e", meanwtim)                << " )  sec" << std::endl
              << "[Min,Max]TimeInMatrixElems  = [ " << string_format("%1.15e", minwtim) << " ,  " << string_format("%1.15e", maxwtim) << " ]  sec" << std::endl
              //<< "StdDevTimeInMatrixElems     = ( " << string_format("%1.15e", stdwtim)                 << " )  sec" << std::endl
              << "TotalTime[MECalcOnly]  (3a) = ( " << string_format("%1.15e", sumw3atim)               << " )  sec" << std::endl
              << "MeanTimeInMECalcOnly        = ( " << string_format("%1.15e", meanw3atim)              << " )  sec" << std::endl
              << "[Min,Max]TimeInMECalcOnly   = [ " << string_format("%1.15e", minw3atim) << " ,  " << string_format("%1.15e", maxw3atim) << " ]  sec" << std::endl
              //<< "StdDevTimeInMECalcOnly      = ( " << string_format("%1.15e", stdw3atim)               << " )  sec" << std::endl
              << std::string(SEP79, '-') << std::endl
              //<< "ProcessID:                  = " << getpid() << std::endl
              //<< "NProcesses                  = " << process.nprocesses << std::endl
              << "TotalEventsComputed         = " << nevtALL << std::endl
              << "EvtsPerSec[Rnd+Rmb+ME](123) = ( " << string_format("%1.15e", nevtALL/(sumgtim+sumrtim+sumwtim)) << " )  sec^-1" << std::endl
              << "EvtsPerSec[Rmb+ME]     (23) = ( " << string_format("%1.15e", nevtALL/(sumrtim+sumwtim))         << " )  sec^-1" << std::endl
              //<< "EvtsPerSec[RndNumbGen]   (1) = ( " << string_format("%1.15e", nevtALL/sumgtim)                  << " )  sec^-1" << std::endl
              //<< "EvtsPerSec[Rambo]        (2) = ( " << string_format("%1.15e", nevtALL/sumrtim)                  << " )  sec^-1" << std::endl
              << "EvtsPerSec[MatrixElems] (3) = ( " << string_format("%1.15e", nevtALL/sumwtim)                   << " )  sec^-1" << std::endl
              << "EvtsPerSec[MECalcOnly] (3a) = ( " << string_format("%1.15e", nevtALL/sumw3atim)                 << " )  sec^-1" << std::endl
              << "EvtsPerSec[Bridge]          = ( " << string_format("%1.15e", bridge_MEs_per_sec)                << " )  sec^-1" << std::endl;
    std::cout << std::string(SEP79, '*') << std::endl
              << "NumMatrixElems(notAbnormal) = " << nevtALL - nabn << std::endl
              << "MeanMatrixElemValue         = ( " << string_format("%1.15e", meanelem) << " +- " << string_format("%1.15e", stdelem/sqrt(nevtALL - nabn)) << " )  GeV^" << meGeVexponent << std::endl // standard error
              << "[Min,Max]MatrixElemValue    = [ " << string_format("%1.15e", minelem)  << " ,  " << string_format("%1.15e", maxelem)                      << " ]  GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue       = ( " << string_format("%1.15e", stdelem)  << " )  GeV^" << meGeVexponent << std::endl
              << "MeanWeight                  = ( " << string_format("%1.15e", meanweig) << " +- " << string_format("%1.15e", stdweig/sqrt(nevtALL - nabn)) << " )" << std::endl // standard error
              << "[Min,Max]Weight             = [ " << string_format("%1.15e", minweig)  << " ,  " << string_format("%1.15e", maxweig)                      << " ]" << std::endl
              << "StdDevWeight                = ( " << string_format("%1.15e", stdweig)  << " )" << std::endl;
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
             << "\"NumIterations\": "      << niter      << ", " << std::endl
             << "\"NumThreadsPerBlock\": " << team_size << ", " << std::endl
             << "\"NumBlocksPerGrid\": " << league_size << ", " << std::endl
             << "\"FP precision\": " << "\"" << FPTYPE_NAME << " (NaN/abnormal=" << nabn << ")\"," << std::endl
             << "\"Complex type\": " << "\"" << COMPLEX_TYPE_NAME << "\"," << std::endl
             //<< "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\"" << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             //<< "\"Momenta memory layout\": " << "\"AOSOA[" << neppM << "]\"" << ( neppM == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Curand generation\": "     << "\"COMMON RANDOM (C++ code)\"," << std::endl;

    jsonFile << "\"NumberOfEntries\": " << niter << "," << std::endl
             << "\"TotalTime[Rnd+Rmb+ME] (123)\": \"" << string_format("%1.15e", sumgtim+sumrtim+sumwtim) << " sec\"," << std::endl
             << "\"TotalTime[Rambo+ME] (23)\": \""    << string_format("%1.15e", sumrtim+sumwtim)         << " sec\"," << std::endl
             << "\"TotalTime[RndNumGen] (1)\": \""    << string_format("%1.15e", sumgtim)                 << " sec\"," << std::endl
             << "\"TotalTime[Rambo] (2)\": \""        << string_format("%1.15e", sumrtim)                 << " sec\"," << std::endl
             << "\"TotalTime[MatrixElems] (3)\": \""  << string_format("%1.15e", sumwtim)                 << " sec\"," << std::endl
             << "\"TotalTime[Bridge]\": \""           << string_format("%1.15e", bridge_time)             << " sec\"," << std::endl
             << "\"MeanTimeInMatrixElems\": \""       << string_format("%1.15e", meanwtim)                << " sec\"," << std::endl
             << "\"MinTimeInMatrixElems\": \""        << std::to_string(minwtim) << " sec\"," << std::endl
             << "\"MaxTimeInMatrixElems\": \""        << std::to_string(maxwtim) << " sec\"," << std::endl
             //<< "ProcessID:                = "        << getpid() << std::endl
             //<< "NProcesses                = "        << process.nprocesses << std::endl
             << "\"TotalEventsComputed\": " << nevtALL << "," << std::endl
             << "\"EvtsPerSec[Rnd+Rmb+ME](123)\": \""  << string_format("%1.15e", nevtALL/(sumgtim+sumrtim+sumwtim)) << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[Rmb+ME] (23)\": \""      << string_format("%1.15e", nevtALL/(sumrtim+sumwtim))         << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[MatrixElems] (3)\": \""  << string_format("%1.15e", nevtALL/sumwtim)                   << " sec^-1\"," << std::endl
             << "\"EvtsPerSec[Bridge]\": \""           << string_format("%1.15e", bridge_MEs_per_sec)                << " sec^-1\"," << std::endl
             << "\"NumMatrixElems(notAbnormal)\": "    << nevtALL - nabn << "," << std::endl
             << "\"MeanMatrixElemValue\": \""          << string_format("%1.15e", meanelem)              << " GeV^" << meGeVexponent << "\"," << std::endl
             << "\"StdErrMatrixElemValue\": \""        << string_format("%1.15e", stdelem/sqrt(nevtALL)) << " GeV^" << meGeVexponent << "\"," << std::endl
             << "\"StdDevMatrixElemValue\": \""        << string_format("%1.15e", stdelem)               << " GeV^" << meGeVexponent << "\"," << std::endl
             << "\"MinMatrixElemValue\": \""           << string_format("%1.15e", minelem)               << " GeV^" << meGeVexponent << "\"," << std::endl
             << "\"MaxMatrixElemValue\": \""           << string_format("%1.15e", maxelem)               << " GeV^" << meGeVexponent << "\"," << std::endl;

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

  } Kokkos::finalize();

  //std::cout << "ALL OK" << std::endl;
  return 0;
}
