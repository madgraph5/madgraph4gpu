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
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "rambo.h"
#include "CommonRandomNumbers.h"

#include "CPPProcess.h"
#include "Memory.h"
#include "timermap.h"

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

bool check_digits(std::string s) {
    return std::all_of(
        s.begin(), s.end(),
        [](char c) { return isdigit(static_cast<unsigned char>(c)); }
    );
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
void debug_me_is_abnormal( const fptype& me, size_t ievtALL )
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
  std::cout << "The '--device_info` flag prints information for all available devices. If a device is chosen by '--device_id', only information for that device is shown." << std::endl;
  std::cout << "The '--device_id' argument selects the device to run code on. (default: 0)" << std::endl;
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
  bool verbose = false;
  bool debug = false;
  bool perf = false;
  bool json = false;
  size_t niter = 0;
  size_t gpublocks = 1;
  size_t gputhreads = 32;
  size_t jsondate = 0;
  size_t jsonrun = 0;
  size_t numvec[5] = {0,0,0,0,0};
  size_t nnum = 0;
  std::string param_card = "../../Cards/param_card.dat";
  bool json_file_bool = false;
  std::string json_file = "";
  std::string d_id_str;
  size_t d_id = 0; // default to device 0
  bool device_chosen = false;
  bool device_info = false;
  auto devices = sycl::device::get_devices();
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
      std::cout << "--bridge unused for now." << std::endl;
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
          for (size_t i=0; i < devices.size(); i++)
          {
            const auto &device = devices[i];
            auto d_name = device.get_info<sycl::info::device::name>();
            std::cout << "    " << i << ": " << d_name << std::endl;
          }
          std::cout << "Terminating. Exit Code: 2" << std::endl;
          return 2;
        }
      }
      else
      {
        std::cout << std::endl << "Invalid device_id, must be integer. Terminating. Exit Code: 3" << std::endl << std::endl;
        return usage(argv[0], 3);
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
      std::cout << "Terminating. Exit Code: 4" << std::endl;
      return 4;
    } else {
      for (size_t i=0; i < devices.size(); i++) {
        std::cout << "device_id " << i << ":" << std::endl;
        print_device_info(devices[i]);
      }
      std::cout << "Terminating. Exit Code: 5" << std::endl;
      return 5;
    }
  }

  static constexpr size_t neppR = mgOnGpu::neppR; // AOSOA layout: constant at compile-time
  if ( gputhreads%neppR != 0 )
  {
    std::cout << "ERROR! #threads/block should be a multiple of neppR=" << neppR << std::endl;
    return usage(argv[0]);
  }

  static constexpr size_t ntpbMAX = mgOnGpu::ntpbMAX;
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

  const size_t ndim = gpublocks * gputhreads; // number of threads in one GPU grid
  const size_t nevt = ndim; // number of events in one iteration == number of GPU threads
  const size_t nevtALL = niter*nevt; // total number of ALL events in all iterations

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

#ifndef MGONGPU_HARDCODE_PARAM
  // Create a process object
  Proc::CPPProcess process( niter, gpublocks, gputhreads, verbose );

  // Read param_card and set parameters
  process.initProc(param_card);
#endif
  const fptype energy = 1500; // historical default, Ecms = 1500 GeV = 1.5 TeV (above the Z peak)
  //const fptype energy = 91.2; // Ecms = 91.2 GeV (Z peak)
  //const fptype energy = 0.100; // Ecms = 100 MeV (well below the Z peak, pure em scattering)
  const int meGeVexponent = -(2 * mgOnGpu::npar - 8);

  // --- 0b. Allocate memory structures
  const std::string alloKey = "0b MemAlloc";
  timermap.start( alloKey );

  // Memory structures for random numbers, momenta, matrix elements and weights on host and device
  static constexpr size_t np4 =  mgOnGpu::np4;
  static constexpr size_t nparf = mgOnGpu::nparf;
  static constexpr size_t npar = mgOnGpu::npar;
  static constexpr size_t ncomb = mgOnGpu::ncomb; // Number of helicity combinations
  const size_t nRnarray = np4*nparf*nevt; // (NB: AOSOA layout with nevt=npagR*neppR events per iteration)
  const size_t nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for AOSOA layouts)
  const size_t nWeights = nevt;
  const size_t nMEs     = nevt; // FIXME: assume process.nprocesses == 1 (eventually: nMEs = nevt * nprocesses?)

  auto hstRnarray   = hstMakeUnique<fptype   >( nRnarray ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  auto hstMomenta   = hstMakeUnique<fptype   >( nMomenta );   // AOSOA[npagM][npar][np4][neppM] (NB: nevt=npagM*neppM)
  auto hstMomenta2  = hstMakeUnique<fptype   >( nMomenta );   // AOSOA[npagM][npar][np4][neppM] (NB: nevt=npagM*neppM)
  auto hstIsGoodHel = hstMakeUnique<bool     >( ncomb );
  auto hstWeights   = hstMakeUnique<fptype   >( nWeights );
  auto hstMEs       = hstMakeUnique<fptype   >( nMEs ); // AOSOA[npagM][neppM] (NB: nevt=npagM*neppM)

  auto devRnarray   = sycl::malloc_device<fptype    >( nRnarray,   q ); 
  auto devMomentaRn = sycl::malloc_device<fptype    >( nMomenta,   q );
  auto devMomenta   = sycl::malloc_device<vector4   >( nMomenta/MGONGPU_VEC_DIM/MGONGPU_FOURVECTOR_DIM,   q );
  auto devRndHel    = sycl::malloc_device<fptype_sv >( nMEs/MGONGPU_VEC_DIM,       q );
  auto devRndCol    = sycl::malloc_device<fptype_sv >( nMEs/MGONGPU_VEC_DIM,       q );
  auto devIsGoodHel = sycl::malloc_device<bool      >( ncomb,      q );
  auto devWeights   = sycl::malloc_device<fptype    >( nWeights,   q );
  auto devMEs       = sycl::malloc_device<fptype_sv >( nMEs/MGONGPU_VEC_DIM,       q );
  auto devSelHel    = sycl::malloc_device<int_sv    >( nMEs/MGONGPU_VEC_DIM,       q );
  auto devSelCol    = sycl::malloc_device<int_sv    >( nMEs/MGONGPU_VEC_DIM,       q );
#ifndef MGONGPU_HARDCODE_PARAM
  auto dev_independent_couplings  = sycl::malloc_device<cxtype>( Proc::independentCouplings::nicoup, q );
  auto dev_independent_parameters = sycl::malloc_device<fptype>( mgOnGpu::nparams,          q );
#endif
  auto devcNGoodHel = sycl::malloc_device<size_t   >( 1,          q ); 
  auto devcGoodHel  = sycl::malloc_device<size_t   >( ncomb,      q ); 

#ifndef MGONGPU_HARDCODE_PARAM
  q.memcpy( dev_independent_couplings, process.get_tIPC_ptr(), 2*Proc::independentCouplings::nicoup*sizeof(fptype) );
  q.memcpy( dev_independent_parameters, process.get_tIPD_ptr(), mgOnGpu::nparams*sizeof(fptype) ).wait();
#endif

  const size_t nbytesRnarray   = nRnarray * sizeof(fptype);
  const size_t nbytesMomenta   = nMomenta * sizeof(fptype);
  const size_t nbytesIsGoodHel = ncomb    * sizeof(bool);
  const size_t nbytesWeights   = nWeights * sizeof(fptype);
  const size_t nbytesMEs       = nMEs     * sizeof(fptype);
  const size_t nbytescNGoodHel =            sizeof(size_t);
  const size_t nbytescGoodHel  = ncomb    * sizeof(size_t);

  constexpr fptype fixedG = 1.2177157847767195; // fixed G for aS=0.118 (hardcoded for now in check_sa.cc, fcheck_sa.f, runTest.cc)
  std::unique_ptr<double[]> genrtimes( new double[niter] );
  std::unique_ptr<double[]> rambtimes( new double[niter] );
  std::unique_ptr<double[]> wavetimes( new double[niter] );
  std::unique_ptr<double[]> wv3atimes( new double[niter] );
  std::unique_ptr<fptype[]> matrixelementALL( new fptype[nevtALL] ); // FIXME: assume process.nprocesses == 1
  std::unique_ptr<fptype[]> weightALL( new fptype[nevtALL] );

  double bridge_time = 0.0;
  #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      mg5amcGpu::Bridge<fptype>* pbridge;
      fptype* _gs;
      fptype* _mes;
      fptype* _rndhel;
      fptype* _rndcol;
      int* _selhel;
      int* _selcol;
      //fptype* _momenta;
      if (bridge) {
          //Initialize bridge object
          pbridge = new mg5amcGpu::Bridge<fptype>( nevt, npar, np4 );
      
          _gs      = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _mes     = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _rndhel  = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _rndcol  = reinterpret_cast<fptype*>(std::malloc(nevt*sizeof(fptype)));
          _selhel  = reinterpret_cast<int*>(std::malloc(nevt*sizeof(int)));
          _selcol  = reinterpret_cast<int*>(std::malloc(nevt*sizeof(int)));
          //_momenta = reinterpret_cast<fptype*>(std::malloc(nevt*npar*np4*sizeof(fptype)));
          //for (size_t i=0; i < nevt*npar*np4; i++) {
          //    _momenta[i] = hstMomenta[i];
          //}
          for (size_t i=0; i < nevt; i++) {
              _gs[i] = fixedG;
          }
      }
  #endif

  // --- 0c. Create curand or common generator
  const std::string cgenKey = "0c GenCreat";
  timermap.start( cgenKey );

  std::vector<std::promise<std::vector<fptype>>> commonRandomPromises;
  CommonRandomNumbers::startGenerateAsync(commonRandomPromises, nRnarray, niter);

  // **************************************
  // *** START MAIN LOOP ON #ITERATIONS ***
  // **************************************

  for (size_t iiter = 0; iiter < niter; ++iiter)
  {
    //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;

    // === STEP 1 OF 3

    // *** START THE OLD-STYLE TIMER FOR RANDOM GEN ***
    double genrtime = 0;

    // --- 1b. Generate all relevant numbers to build nevt events (i.e. nevt phase space points) on the host
    const std::string rngnKey = "1b GenRnGen";
    timermap.start( rngnKey );

    std::vector<fptype> commonRnd = commonRandomPromises[iiter].get_future().get();
    assert( nRnarray == static_cast<size_t>( commonRnd.size() ) );

    memcpy( hstRnarray.get(), commonRnd.data(), nRnarray * sizeof(fptype) );

    // --- 1c. Copy rnarray from host to device
    const std::string htodKey = "1c CpHTDrnd";
    genrtime += timermap.start( htodKey );

    q.memcpy(devRnarray, hstRnarray.get(), nbytesRnarray).wait();

    // *** STOP THE OLD-STYLE TIMER FOR RANDOM GEN ***
    genrtime += timermap.stop();

    // === STEP 2 OF 3
    // Fill in particle momenta for each of nevt events on the device

    // *** START THE OLD-STYLE TIMER FOR RAMBO ***
    double rambtime = 0;

    // --- 2a. Fill in momenta of initial state particles on the device
    const std::string riniKey = "2a RamboIni";
    timermap.start( riniKey );

    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{gpublocks}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                rambo2toNm0::ramboGetMomentaInitial( energy, devMomentaRn, ievt );
            });
        }));
    });
    q.wait();

    // --- 2b. Fill in momenta of final state particles using the RAMBO algorithm on the device
    // (i.e. map random numbers to final-state particle momenta for each of nevt events)
    const std::string rfinKey = "2b RamboFin";
    rambtime += timermap.start( rfinKey );

    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for_work_group(sycl::range<1>{gpublocks}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);
                rambo2toNm0::ramboGetMomentaFinal( energy, devRnarray, devMomentaRn, devWeights, ievt );
            });
        }));
    });
    q.wait();

    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    rambtime += timermap.start( cwgtKey );
    q.memcpy(hstWeights.get(), devWeights, nbytesWeights).wait();

    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    rambtime += timermap.start( cmomKey );
    q.memcpy(hstMomenta.get(), devMomentaRn, nbytesMomenta).wait();
    mg5amcGpu::hst_transposeMomenta(hstMomenta2.get(), hstMomenta.get(), nevt);
    q.memcpy(devMomenta, hstMomenta2.get(), nbytesMomenta).wait();

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
      q.submit([&](sycl::handler& cgh) {
          #if MGONGPU_VEC_DIM > 1
          cgh.parallel_for_work_group(sycl::range<1>{gpublocks/MGONGPU_VEC_DIM}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
          #else
          cgh.parallel_for_work_group(sycl::range<1>{gpublocks}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
          #endif
              #ifndef MGONGPU_HARDCODE_PARAM
                  //Load independent couplings and parameters into shared memory if not hardcoded
                  auto _dev_independent_couplings = dev_independent_couplings;
                  auto dev_parameters = dev_independent_parameters;
              #endif
              wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                  size_t ievt = index.get_global_id(0);

                  #ifdef MGONGPU_HARDCODE_PARAM
                      //Load parameters into local (private) memory if hardcoded
                      auto dev_parameters = Proc::independent_parameters<fptype>;
                  #endif

                  //Load helicities and couplings into local (private) memory
                  auto dev_helicities = Proc::helicities<signed char>;
                  cxtype_sv dev_couplings[Proc::dependentCouplings::ndcoup + Proc::independentCouplings::nicoup];

                  #if MGONGPU_NDCOUP > 0
                      Proc::dependentCouplings::set_couplings_from_G(dev_couplings, fptype_sv(fixedG)); 
                  #endif

                  #if MGONGPU_NICOUP > 0
                      #ifdef MGONGPU_HARDCODE_PARAM
                          //Load independent couplings into local (private) memory if hardcoded
                          auto _dev_independent_couplings = Proc::independentCouplings::independent_couplings<cxtype, fptype>;
                      #endif
                      for (size_t i = 0; i < Proc::independentCouplings::nicoup; i++) {
                          dev_couplings[Proc::dependentCouplings::ndcoup + i] = _dev_independent_couplings[i];
                      }
                  #endif

                  Proc::sigmaKin_getGoodHel( devMomenta + npar*ievt, devIsGoodHel, dev_helicities, dev_couplings, dev_parameters );
              });
          }));
      });
      q.wait();

      // ... 0d2. Copy back good helicity mask to the host
      q.memcpy(hstIsGoodHel.get(), devIsGoodHel, nbytesIsGoodHel).wait();

      // ... 0d3. Copy back good helicity list to constant memory on the device
      size_t goodHel[mgOnGpu::ncomb] = {0};
      size_t nGoodHel = Proc::sigmaKin_setGoodHel( hstIsGoodHel.get(), goodHel );

      q.memcpy( devcNGoodHel, &nGoodHel, nbytescNGoodHel ).wait();
      q.memcpy( devcGoodHel, goodHel, nbytescGoodHel ).wait();
    }

    // *** START THE OLD-STYLE TIMERS FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    double wavetime = 0; // calc plus copy
    double wv3atime = 0; // calc only

    // --- 3a. SigmaKin
    const std::string skinKey = "3a SigmaKin";
    timermap.start( skinKey );

    q.submit([&](sycl::handler& cgh) {
        #if MGONGPU_VEC_DIM > 1
        cgh.parallel_for_work_group(sycl::range<1>{gpublocks/MGONGPU_VEC_DIM}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
        #else
        cgh.parallel_for_work_group(sycl::range<1>{gpublocks}, sycl::range<1>{gputhreads}, ([=](sycl::group<1> wGroup) {
        #endif
            #ifndef MGONGPU_HARDCODE_PARAM
                //Load independent couplings and parameters into shared memory if not hardcoded
                auto _dev_independent_couplings = dev_independent_couplings;
                auto dev_parameters = dev_independent_parameters;
            #endif
            wGroup.parallel_for_work_item([&](sycl::h_item<1> index) {
                size_t ievt = index.get_global_id(0);

                #ifdef MGONGPU_HARDCODE_PARAM
                    //Load parameters into local (private) memory if hardcoded
                    auto dev_parameters = Proc::independent_parameters<fptype>;
                #endif

                //Load helicities and couplings into local (private) memory
                auto dev_helicities = Proc::helicities<signed char>;
                cxtype_sv dev_couplings[Proc::dependentCouplings::ndcoup + Proc::independentCouplings::nicoup];

                #if MGONGPU_NDCOUP > 0
                    Proc::dependentCouplings::set_couplings_from_G(dev_couplings, fptype_sv(fixedG)); 
                #endif

                #if MGONGPU_NICOUP > 0
                    #ifdef MGONGPU_HARDCODE_PARAM
                        //Load independent couplings into local (private) memory if hardcoded
                        auto _dev_independent_couplings = Proc::independentCouplings::independent_couplings<cxtype, fptype>;
                    #endif
                    for (size_t i = 0; i < Proc::independentCouplings::nicoup; i++) {
                        dev_couplings[Proc::dependentCouplings::ndcoup + i] = _dev_independent_couplings[i];
                    }
                #endif

                #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                    devMEs[ievt] = Proc::sigmaKin( devMomenta + npar*ievt, devRndHel + ievt, devRndCol + ievt, devSelHel + ievt, devSelCol + ievt, 0, dev_helicities, dev_couplings, dev_parameters, devcNGoodHel, devcGoodHel );
                #else
                    devMEs[ievt] = Proc::sigmaKin( devMomenta + npar*ievt, devRndHel + ievt, devRndCol + ievt, devSelHel + ievt, devSelCol + ievt, dev_helicities, dev_couplings, dev_parameters, devcNGoodHel, devcGoodHel );
                #endif
            });
        }));
    });
    q.wait();

    // *** STOP THE NEW OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wv3atime += timermap.stop(); // calc only
    wavetime += wv3atime; // calc plus copy

    // --- 3b. CopyDToH MEs
    const std::string cmesKey = "3b CpDTHmes";
    timermap.start( cmesKey );

    q.memcpy(hstMEs.get(), devMEs, nbytesMEs).wait();

    // *** STOP THE OLD OLD-STYLE TIMER FOR MATRIX ELEMENTS (WAVEFUNCTIONS) ***
    wavetime += timermap.stop(); // calc plus copy

    #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        if (bridge) {
            const std::string bridgeKey = "3c bridge";
            timermap.start( bridgeKey );
            pbridge->gpu_sequence( hstMomenta.get(), _gs, _rndhel, _rndcol, 0, _mes, _selhel, _selcol );
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
                    << std::setw(14) << hstMomenta[npar*np4*ievt + ipar*np4 + 0] // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << hstMomenta[npar*np4*ievt + ipar*np4 + 1] // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << hstMomenta[npar*np4*ievt + ipar*np4 + 2] // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << hstMomenta[npar*np4*ievt + ipar*np4 + 3] // AOSOA[ipagM][ipar][3][ieppM]
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

  }

  #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if (bridge) {
          // Destroy bridge objects
          std::free(_gs);
          std::free(_mes);
          std::free(_rndhel);
          std::free(_rndcol);
          std::free(_selhel);
          std::free(_selcol);
          //std::free(_momenta);
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
  //double sqsgtim = 0;
  double mingtim = genrtimes[0];
  double maxgtim = genrtimes[0];
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
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
  for ( size_t iiter = 0; iiter < niter; ++iiter )
  {
    sumw3atim += wv3atimes[iiter];
    //sqsw3atim += wv3atimes[iiter]*wv3atimes[iiter];
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
  {
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
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
  for ( size_t ievtALL = 0; ievtALL < nevtALL; ++ievtALL )
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
  const std::string syclfrKey = "9a sycl_free";
  timermap.start( syclfrKey );

  sycl::free( devRnarray   , q);
  sycl::free( devMomenta   , q);
  sycl::free( devIsGoodHel , q);
  sycl::free( devWeights   , q);
  sycl::free( devMEs       , q);
  sycl::free( devcNGoodHel , q);
  sycl::free( devcGoodHel  , q);

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
              << "Process                     = " << XSTRINGIFY(MG_EPOCH_PROCESS_ID) << "_SYCL"
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
              << "Complex type                = EXTRA" << std::endl
              << "RanNumb memory layout       = AOSOA[" << neppR << "]" << ( neppR == 1 ? " == AOS" : "" ) << " [HARDCODED FOR REPRODUCIBILITY]" << std::endl
              << "Random number generation    = COMMON RANDOM (C++ code)" << std::endl
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
             << "\"NumThreadsPerBlock\": " << gputhreads << ", " << std::endl
             << "\"NumBlocksPerGrid\": "   << gpublocks  << ", " << std::endl
#if defined MGONGPU_FPTYPE_DOUBLE
             << "\"FP precision\": " << "\"DOUBLE (NaN/abnormal=" << nabn << ")\"," << std::endl
#elif defined MGONGPU_FPTYPE_FLOAT
             << "\"FP precision\": " << "\"FLOAT (NaN/abnormal=" << nabn << ")\"," << std::endl
#endif
             << "\"Complex type\": " << "\"EXTRA\"," << std::endl
             << "\"RanNumb memory layout\": " << "\"AOSOA[" << neppR << "]\"" << ( neppR == 1 ? " == AOS" : "" ) << ", " << std::endl
             << "\"Curand generation\": " << "\"COMMON RANDOM (C++ code)\"," << std::endl;

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

  //std::cout << "ALL OK" << std::endl;
  return 0;
}
