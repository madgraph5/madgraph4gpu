#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <unistd.h>

#include <CL/sycl.hpp>
#include "mgOnGpuConfig.h"

#include "CPPProcess.h"
#include "CrossSectionKernels.h"
#include "MatrixElementKernels.h"
#include "MemoryBuffers.h"
#include "RamboSamplingKernels.h"
#include "RandomNumberKernels.h"
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

bool check_digits(std::string s) {
  return std::all_of(
    s.begin(), s.end(),
    [](char c) { return isdigit(static_cast<unsigned char>(c)); }
  );
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j] [--param_card <PARAM_CARD_FILE>] [--json_file <JSON_FILE>] [--bridge] [--device_info] --device_id <DEVICE_ID>"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only enables NaN/abnormal warnings and OMP debugging" << std::endl;
  std::cout << "The '--device_info` flag prints information for all available devices. If a device is chosen by '--device_id', only information for that device is shown." << std::endl;
  std::cout << "The '--device_id' arguments selects the device to run code on. (default: 0)" << std::endl;
  std::cout << "The '--help|-h' flag prints this message" << std::endl;
  return ret;
}

void print_device_type( const sycl::info::device_type dt ) {
    if (dt == sycl::info::device_type::cpu) {
        std::cout << "cpu"; }
    if (dt == sycl::info::device_type::gpu) {
        std::cout << "gpu"; }
    if (dt == sycl::info::device_type::accelerator) {
        std::cout << "accelerator"; }
    if (dt == sycl::info::device_type::custom) {
        std::cout << "custom"; }
    if (dt == sycl::info::device_type::automatic) {
        std::cout << "automatic"; }
    if (dt == sycl::info::device_type::host) {
        std::cout << "host"; }
    if (dt == sycl::info::device_type::all) {
        std::cout << "all"; }
}

void print_device_info(const sycl::device& device) {
    std::cout << "    name: " << device.get_info<sycl::info::device::name>() << std::endl;
    std::cout << "    platform: " << device.get_info<sycl::info::device::platform>().get_info<sycl::info::platform::name>() << std::endl;
    std::cout << "    vendor: " << device.get_info<sycl::info::device::vendor>() << std::endl;
    std::cout << "    vendor_id: " << device.get_info<sycl::info::device::vendor_id>() << std::endl;
    std::cout << "    driver_version: " << device.get_info<sycl::info::device::driver_version>() << std::endl;
    std::cout << "    global_mem_size: " << device.get_info<sycl::info::device::global_mem_size>() << std::endl;
    std::cout << "    local_mem_size: " << device.get_info<sycl::info::device::local_mem_size>() << std::endl;

    auto workgroup_size = device.get_info<sycl::info::device::max_work_group_size>();
    auto max_compute_units = device.get_info<sycl::info::device::max_compute_units>();
    //auto n_groups = (num_steps - 1) / workgroup_size + 1;
    //n_groups = std::min(decltype(n_groups)(max_compute_units),n_groups);  // make groups max number of compute units or less
    std::cout << "    workgroup_size: " << workgroup_size << std::endl;
    std::cout << "    max_compute_units: " << max_compute_units << std::endl;

    std::cout << "    usm_support: ";
    if (device.has(sycl::aspect::usm_device_allocations)) {
        std::cout << "yes" << std::endl;
    } else {
        std::cout << "no" << std::endl;
    }
    std::cout << "    device_type: ";
    print_device_type(device.get_info<sycl::info::device::device_type>());
    std::cout << std::endl;
}

void print_device_info() {
    auto platforms = sycl::platform::get_platforms();
    for (const auto &platform: platforms) {
        std::cout << platform.get_backend() << ':' << std::endl;
        std::cout << platform.get_info<sycl::info::platform::name>() << ':' << std::endl;
        for (const auto &device: platform.get_devices()) {
            print_device_info(device);
        }
    }
    std::cout << std::endl;
}

int main(int argc, char **argv)
{
  // Namespaces for CUDA and C++ (FIXME - eventually use the same namespace everywhere...)
  using namespace mg5amcGpu;
  
  // DEFAULTS FOR COMMAND LINE ARGUMENTS
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
  std::string d_id_str;
  size_t d_id = 0; // default to device 0
  bool device_chosen = false;
  bool device_info = false;
  auto devices = sycl::device::get_devices();
  // Random number mode
  enum class RandomNumberMode{ CommonRandom=0, CurandHost=1, CurandDevice=2 };
  RandomNumberMode rndgen = RandomNumberMode::CommonRandom; // default on CPU if build has no curand
  // Rambo sampling mode (NB RamboHost implies CommonRandom or CurandHost!)
  enum class RamboSamplingMode{ RamboHost=1, RamboDevice=2 };
  RamboSamplingMode rmbsmp = RamboSamplingMode::RamboHost; // default on CPU
  // Bridge emulation mode (NB Bridge implies RamboHost!)
  bool bridge = false;

  // READ COMMAND LINE ARGUMENTS
  for ( int argn = 1; argn < argc; ++argn )
  {
    std::string arg = argv[argn];
    if ( ( arg == "--verbose" ) || ( arg == "-v" ) )
    {
      verbose = true;
    }
    else if ( ( arg == "--debug" ) || ( arg == "-d" ) )
    {
      debug = true;
    }
    else if ( ( arg == "--performance" ) || ( arg == "-p" ) )
    {
      perf = true;
    }
    else if ( ( arg == "--json" ) || ( arg == "-j" ) )
    {
      json = true;
    }
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
    else if ( arg == "--device_info" )
    {
        device_info = true;
    }
    else if ( arg == "--device_id" )
    {
      d_id_str = argv[argn + 1];
      if (check_digits(d_id_str))
      {
        d_id = std::stoi(d_id_str);
        device_chosen = true;
        if (d_id >= devices.size())
        {
          std::cout << "Invalid device_id. Please choose device from: " << std::endl;
          for (int i=0; i < devices.size(); i++)
          {
            const auto &device = devices[i];
            auto d_name = device.get_info<sycl::info::device::name>();
            std::cout << "    " << i << ": " << d_name << std::endl;
          }
          std::cout << "Terminating. Exit Code: -999" << std::endl;
          return -999;
        }
      }
      else
      {
        std::cout << std::endl << "Invalid device_id, must be integer. Terminating. Exit Code: -998" << std::endl << std::endl;
        return usage(argv[0], -998);
      }
      argn++;
    }
    else if ( is_number(argv[argn]) && nnum<5 )
    {
      numvec[nnum++] = atoi( argv[argn] );
    }
    else
    {
      return usage( argv[0] );
    }
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

  if (device_info) {
    if (device_chosen) {
      print_device_info(devices[d_id]);
      std::cout << "Terminating. Exit Code: -997" << std::endl;
      return -997;
    } else {
      for (int i=0; i < devices.size(); i++) {
        std::cout << "device_id " << i << ":" << std::endl;
        print_device_info(devices[i]);
      }
      std::cout << "Terminating. Exit Code: -996" << std::endl;
      return -996;
    }
  }

  constexpr int neppM = mgOnGpu::neppM; // AOSOA layout
  constexpr int neppR = 8; // FIXME HARDCODED

  using mgOnGpu::ntpbMAX;
  if ( gputhreads > ntpbMAX )
  {
    std::cout << "ERROR! #threads/block should be <= " << ntpbMAX << std::endl;
    return usage(argv[0]);
  }

  const int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  const int nevt = ndim; // number of events in one iteration == number of GPU threads

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // *** START THE NEW TIMERS ***
  mgOnGpu::TimerMap timermap;

  // === STEP 0 - INITIALISE

  // --- 00. Initialise SYCL
  sycl::queue q = sycl::queue(devices[d_id]);

  if (verbose) {
    auto device = q.get_device();
    std::cout << "Selected " << device.get_info<sycl::info::device::name>()
              << " on platform "
              << device.get_info<sycl::info::device::platform>().get_info<sycl::info::platform::name>()
              << std::endl;
  }

  // --- 0a. Initialise physics process
  const std::string procKey = "0a ProcInit";
  timermap.start( procKey );

  // Create a process object
  CPPProcess process( niter, gpublocks, gputhreads, verbose );

  // Read param_card and set parameters
  process.initProc(param_card);
  const fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
  //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
  const int meGeVexponent = -(2 * mgOnGpu::npar - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory buffers for random numbers
  PinnedHostBufferRandomNumbers hstRnarray( nevt, q );
  DeviceBufferRandomNumbers devRnarray( nevt, q );

  // Memory buffers for momenta
  PinnedHostBufferMomenta hstMomenta( nevt, q );
  DeviceBufferMomenta devMomenta( nevt, q );

  // Memory buffers for sampling weights
  PinnedHostBufferWeights hstWeights( nevt, q );
  DeviceBufferWeights devWeights( nevt, q );

  // Memory buffers for matrix elements
  PinnedHostBufferMatrixElements hstMatrixElements( nevt, q );
  DeviceBufferMatrixElements devMatrixElements( nevt, q );

  // Memory buffers for the helicity mask
  using mgOnGpu::ncomb; // the number of helicity combinations
  //FIXME SYCL need queue for these variables
  PinnedHostBufferHelicityMask hstIsGoodHel( ncomb, q );
  DeviceBufferHelicityMask devIsGoodHel( ncomb, q );

  std::unique_ptr<double[]> genrtimes( new double[niter] );
  std::unique_ptr<double[]> rambtimes( new double[niter] );
  std::unique_ptr<double[]> wavetimes( new double[niter] );
  std::unique_ptr<double[]> wv3atimes( new double[niter] );

  // --- 0c. Create curand or common generator
  const std::string cgenKey = "0c GenCreat";
  timermap.start( cgenKey );
  // Allocate the appropriate RandomNumberKernel
  std::unique_ptr<RandomNumberKernelBase> prnk;
  if ( rndgen == RandomNumberMode::CommonRandom )
  {
    prnk.reset( new CommonRandomNumberKernel( hstRnarray ) );
  }

  // --- 0c. Create rambo sampling kernel [keep this in 0c for the moment]
  std::unique_ptr<SamplingKernelBase> prsk;
  if ( rmbsmp == RamboSamplingMode::RamboHost )
  {
    prsk.reset( new RamboSamplingKernelHost( energy, hstRnarray, hstMomenta, hstWeights, nevt ) );
  }

  // --- 0c. Create matrix element kernel [keep this in 0c for the moment]
  std::unique_ptr<MatrixElementKernelBase> pmek;
  pmek.reset( new MatrixElementKernelDevice( devMomenta, devMatrixElements, q, gpublocks, gputhreads ) );
  pmek->setDeviceArrays( process.get_tHel_ptr(), process.get_tIPC_ptr(), process.get_tIPD_ptr() );

  // --- 0c. Create cross section kernel [keep this in 0c for the moment]
  EventStatistics hstStats;
  CrossSectionKernelHost xsk( hstWeights, hstMatrixElements, hstStats, nevt );

  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************

  for (int iiter = 0; iiter < niter; ++iiter)
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
    const unsigned long long seed = 20200805;
    const std::string sgenKey = "1a GenSeed ";
    timermap.start( sgenKey );
    prnk->seedGenerator( seed+iiter );
    genrtime += timermap.stop();

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );
    prnk->generateRnarray();
    //std::cout << "Got random numbers" << std::endl;

      // --- 1c. Copy rnarray from host to device
      const std::string htodKey = "1c CpHTDrnd";
      genrtime += timermap.start( htodKey );
      copyDeviceFromHost( devRnarray, hstRnarray );

    // *** STOP THE OLD-STYLE TIMER FOR RANDOM GEN ***
    genrtime += timermap.stop();

    // === STEP 2 OF 3
    // Fill in particle momenta for each of nevt events on the device

    // *** START THE OLD-STYLE TIMER FOR RAMBO ***
    double rambtime = 0;

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );
    prsk->getMomentaInitial();
    //std::cout << "Got initial momenta" << std::endl;

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );
    prsk->getMomentaFinal();
    //std::cout << "Got final momenta" << std::endl;

    // --- 2c. CopyHToD Weights
    const std::string cwgtKey = "2c CpHTDwgt";
    rambtime += timermap.start( cwgtKey );
    copyDeviceFromHost( devWeights, hstWeights );
    
    // --- 2d. CopyHToD Momenta
    const std::string cmomKey = "2d CpHTDmom";
    rambtime += timermap.start( cmomKey );
    copyDeviceFromHost( devMomenta, hstMomenta );

    // *** STOP THE OLD-STYLE TIMER FOR RAMBO ***
    rambtime += timermap.stop();

    // === STEP 3 OF 3
    // Evaluate matrix elements for all nevt events
    // 0d. For Bridge only, transpose C2F [renamed as 0d: this is not initialisation, but I want it out of the ME timers (#371)]
    // 0e. (Only on the first iteration) Get good helicities [renamed as 0e: this IS initialisation!]
    // 3a. Evaluate MEs on the device (include transpose F2C for Bridge)
    // 3b. Copy MEs back from device to host

    // --- 0d. TransC2F
    if ( bridge )
    {
      const std::string tc2fKey = "0d TransC2F";
      timermap.start( tc2fKey );
      //dynamic_cast<BridgeKernelBase*>( pmek.get() )->transposeInputMomentaC2F();
    }    

    // --- 0e. SGoodHel
    if ( iiter == 0 )
    {
      const std::string ghelKey = "0e SGoodHel";
      timermap.start( ghelKey );
      pmek->computeGoodHelicities();
    }

    // *** START THE OLD-STYLE TIMERS FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    double wavetime = 0; // calc plus copy
    double wv3atime = 0; // calc only

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );
    pmek->computeMatrixElements();

    // *** STOP THE NEW OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wv3atime += timermap.stop(); // calc only
    wavetime += wv3atime; // calc plus copy

    if ( ! bridge )
    {
      // --- 3b. CopyDToH MEs
      const std::string cmesKey = "3b CpDTHmes";
      timermap.start( cmesKey );
      copyHostFromDevice( hstMatrixElements, devMatrixElements );
      // *** STOP THE OLD OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
      wavetime += timermap.stop(); // calc plus copy
    }

    // === STEP 4 FINALISE LOOP
    // --- 4@ Update event statistics
    const std::string updtKey = "4@ UpdtStat";
    timermap.start(updtKey);
    xsk.updateEventStatistics();

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
        int ipagM = ievt/neppM;
        int ieppM = ievt%neppM;
        // Display momenta
        std::cout << "Momenta:" << std::endl;
        for (int ipar = 0; ipar < mgOnGpu::npar; ipar++)
        {
          // NB: 'setw' affects only the next field (of any type)
          std::cout << std::scientific // fixed format: affects all floats (default precision: 6)
                    << std::setw(4) << ipar + 1
                    << std::setw(14) << hstMomenta.data()[ipagM*mgOnGpu::npar*mgOnGpu::np4*neppM + ipar*mgOnGpu::np4*neppM + 0*neppM + ieppM]
                    << std::setw(14) << hstMomenta.data()[ipagM*mgOnGpu::npar*mgOnGpu::np4*neppM + ipar*mgOnGpu::np4*neppM + 1*neppM + ieppM]
                    << std::setw(14) << hstMomenta.data()[ipagM*mgOnGpu::npar*mgOnGpu::np4*neppM + ipar*mgOnGpu::np4*neppM + 2*neppM + ieppM]
                    << std::setw(14) << hstMomenta.data()[ipagM*mgOnGpu::npar*mgOnGpu::np4*neppM + ipar*mgOnGpu::np4*neppM + 3*neppM + ieppM]
                    << std::endl
                    << std::defaultfloat; // default format: affects all floats
        }
        std::cout << std::string(SEP79, '-') << std::endl;
        // Display matrix elements
        std::cout << " Matrix element = " << hstMatrixElements.data()[ievt]
                  << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(SEP79, '-') << std::endl;
      }
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
  //double sqsgtim = 0;
  double mingtim = genrtimes[0];
  double maxgtim = genrtimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumgtim += genrtimes[iiter];
    //sqsgtim += genrtimes[iiter]*genrtimes[iiter];
    mingtim = std::min( mingtim, genrtimes[iiter] );
    maxgtim = std::max( maxgtim, genrtimes[iiter] );
  }

  double sumrtim = 0;
  //double sqsrtim = 0;
  double minrtim = rambtimes[0];
  double maxrtim = rambtimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumrtim += rambtimes[iiter];
    //sqsrtim += rambtimes[iiter]*rambtimes[iiter];
    minrtim = std::min( minrtim, rambtimes[iiter] );
    maxrtim = std::max( maxrtim, rambtimes[iiter] );
  }

  double sumwtim = 0;
  //double sqswtim = 0;
  double minwtim = wavetimes[0];
  double maxwtim = wavetimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumwtim += wavetimes[iiter];
    //sqswtim += wavetimes[iiter]*wavetimes[iiter];
    minwtim = std::min( minwtim, wavetimes[iiter] );
    maxwtim = std::max( maxwtim, wavetimes[iiter] );
  }
  double meanwtim = sumwtim / niter;
  //double stdwtim = std::sqrt( sqswtim / niter - meanwtim * meanwtim );

  double sumw3atim = 0;
  //double sqsw3atim = 0;
  double minw3atim = wv3atimes[0];
  double maxw3atim = wv3atimes[0];
  for ( int iiter = 0; iiter < niter; ++iiter )
  {
    sumw3atim += wv3atimes[iiter];
    //sqsw3atim += wv3atimes[iiter]*wv3atimes[iiter];
    minw3atim = std::min( minw3atim, wv3atimes[iiter] );
    maxw3atim = std::max( maxw3atim, wv3atimes[iiter] );
  }
  double meanw3atim = sumw3atim / niter;
  //double stdw3atim = std::sqrt( sqsw3atim / niter - meanw3atim * meanw3atim );

  const int nevtALL = hstStats.nevtALL; // total number of ALL events in all iterations
  if ( nevtALL !=  niter*nevt )
    std::cout << "ERROR! nevtALL mismatch " << nevtALL << " != " << niter*nevt << std::endl; // SANITY CHECK
  int nabn = hstStats.nevtABN;
  int nzero = hstStats.nevtZERO;

  // === STEP 9 FINALISE
  
  std::string rndgentxt;
  if ( rndgen == RandomNumberMode::CommonRandom ) rndgentxt = "COMMON RANDOM HOST";

  // Workflow description summary
  std::string wrkflwtxt;
  wrkflwtxt += "SYC:";
  // -- DOUBLE or FLOAT?
#if defined MGONGPU_FPTYPE_DOUBLE
  wrkflwtxt += "DBL+";
#elif defined MGONGPU_FPTYPE_FLOAT
  wrkflwtxt += "FLT+";
#else
  wrkflwtxt += "???+"; // no path to this statement
#endif
  // -- extra.h complex numbers
  wrkflwtxt += "XTR:";

  // -- COMMON random numbers
  if ( rndgen == RandomNumberMode::CommonRandom ) wrkflwtxt += "COMMON+";
  else wrkflwtxt += "??????+"; // no path to this statement

  // -- HOST rambo sampling
  if ( rmbsmp == RamboSamplingMode::RamboHost ) wrkflwtxt += "RMBHST+";
  else wrkflwtxt += "??????+"; // no path to this statement

  // -- DEVICE matrix elements? Standalone MEs or BRIDGE?
  if ( !bridge ) wrkflwtxt += "MESDEV";
  else wrkflwtxt += "BRDDEV";

  // -- SIMD matrix elements?
  wrkflwtxt += "/none";

  // -- Has cxtype_v::operator[] bracket with non-const reference?
  wrkflwtxt += "+NAVBRK"; // N/A

  // --- 9a Dump to screen
  const std::string dumpKey = "9a DumpScrn";
  timermap.start(dumpKey);

  if (!(verbose || debug || perf))
  {
    std::cout << std::endl;
  }

  if (perf)
  {
  // Dump all configuration parameters and all results
    std::cout << std::string(SEP79, '*') << std::endl
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_SYCL"
              << " [CXX]"
              //<< " [" << process.getCompiler() << "]"
#ifdef MGONGPU_INLINE_HELAMPS
              << " [inlineHel=1]"
#else
              << " [inlineHel=0]"
#endif
              << "NumBlocksPerGrid            = " << gpublocks << std::endl
              << "NumThreadsPerBlock          = " << gputhreads << std::endl
              << "NumIterations               = " << niter << std::endl
              << std::string(SEP79, '-') << std::endl;
    std::cout << "Workflow summary            = " << wrkflwtxt << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
              << "FP precision                = DOUBLE (NaN/abnormal=" << nabn << ", zero=" << nzero << ")" << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
              << "FP precision                = FLOAT (NaN/abnormal=" << nabn << ", zero=" << nzero << ")" << std::endl
#endif
              << "Complex type                = EXTRA::COMPLEX" << std::endl
              << "RanNumb memory layout       = AOSOA[" << neppR << "]"
              << ( neppR == 1 ? " == AOS" : "" )
              << " [HARDCODED FOR REPRODUCIBILITY]" << std::endl
              << "Momenta memory layout       = AOSOA[" << neppM << "]"
              << ( neppM == 1 ? " == AOS" : "" ) << std::endl
              << "Random number generation    = " << rndgentxt << std::endl
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
              << hstStats;
  }

  // --- 9b Dump to json
  const std::string jsonKey = "9b DumpJson";
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
             << "\"EXTRA::COMPLEX\"," << std::endl
             << "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\""
             << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Momenta memory layout\": " << "\"AOSOA[" << neppM << "]\""
             << ( neppM == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Curand generation\": "
             << "\"" << rndgentxt << "\"," << std::endl;

    double minelem = hstStats.minME;
    double maxelem = hstStats.maxME;
    double meanelem = hstStats.meanME();
    double stdelem = hstStats.stdME();

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

  // [NB some resources like curand generators will be deleted here when stack-allocated classes go out of scope]
  //std::cout << "ALL OK" << std::endl;
  return 0;
}
