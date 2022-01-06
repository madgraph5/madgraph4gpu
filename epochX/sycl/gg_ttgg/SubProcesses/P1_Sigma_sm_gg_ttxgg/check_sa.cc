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

#ifdef SYCL_LANGUAGE_VERSION
#include "rambo.cc"
#include "CPPProcess.cc"
#else
#include "rambo.h"
#include "CPPProcess.h"
#endif

#ifdef MGONGPU_COMMONRAND_ONHOST
#include "CommonRandomNumbers.h"
#endif
#include "timermap.h"

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
            << " [--verbose|-v] [--debug|-d] [--performance|-p] [--json|-j] [--param_card <PARAM_CARD_FILE>] [--json_file <JSON_FILE>] --device_id <DEVICE_ID>"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl << std::endl;
  std::cout << "The number of events per iteration is #gpuBlocksPerGrid * #gpuThreadsPerBlock" << std::endl;
  std::cout << "(also in CPU/C++ code, where only the product of these two parameters counts)" << std::endl << std::endl;
  std::cout << "Summary stats are always computed: '-p' and '-j' only control their printout" << std::endl;
  std::cout << "The '-d' flag only controls if nan's emit warnings" << std::endl;
  std::cout << "The '--device_info` flag prints information for all available devices. If a device is chosen by '--device_id', only information for that device is shown." << std::endl;
  std::cout << "The '--device_id' arguments selects the device to run code on. (default: 0)" << std::endl;
  std::cout << "The '--help|-h' flag prints this message" << std::endl;
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
std::unique_ptr<T[]> hstMakeUnique(std::size_t N) { return std::unique_ptr<T[]>{ new T[N] }; };
#endif

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

class user_input_selector : public sycl::device_selector {
    public:
        std::vector<uint32_t> d_ids;
        std::vector<sycl::device> devices;
        uint32_t d_id;
        user_input_selector();
        user_input_selector(uint32_t);
        int operator()(const sycl::device& dev) const override {
            auto device_id = dev.get_info<sycl::info::device::vendor_id>();
            if (device_id == d_id) {
                return 1;
            }
            return -1;
        }
};

user_input_selector::user_input_selector(void) {
    bool device_chosen = false;
    devices = sycl::device::get_devices();
    for (const auto &device: devices) {
        d_ids.push_back(device.get_info<sycl::info::device::vendor_id>());
    }
    print_device_info();
    std::string d_id_str;
    std::cout << "Choose device by entering vendor_id: ";
    std::cin >> d_id_str;
    bool are_digits = check_digits(d_id_str);
    if (are_digits) {
        d_id = std::stoi(d_id_str);
        if (std::find(d_ids.begin(), d_ids.end(), d_id) != d_ids.end()) {
            device_chosen = true;
        }
    }
    while (!device_chosen) {
        std::cout << "Invalid vendor_id. Please choose from ( ";
        for (auto _d_id: d_ids) {
            std::cout << _d_id << " ";
        }
        std::cout << ")." << std::endl;
        std::cout << "Choose device by entering vendor_id: ";
        std::cin >> d_id_str;
        are_digits = check_digits(d_id_str);
        if (are_digits) {
            d_id = std::stoi(d_id_str);
            if (std::find(d_ids.begin(), d_ids.end(), d_id) != d_ids.end()) {
                device_chosen = true;
            }
        }
    }
}

user_input_selector::user_input_selector(uint32_t _d_id) {
    d_id = _d_id;
}

sycl::queue get_sycl_queue(bool device_chosen,uint32_t d_id) {
    if (device_chosen) {
        sycl::queue q( user_input_selector{d_id} );
        return q;
    } else {
        sycl::queue q( user_input_selector{} );
        return q;
    }
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
  std::string d_id_str;
  uint32_t d_id = 0; // default to device 0
  bool device_chosen = false;
  bool device_info = false;
  auto devices = sycl::device::get_devices();
  if (devices.size() == 0) {
      std::cout << "No SYCL devices detected." << std::endl;
      std::cout << "Terminating. Exit Code: -995" << std::endl;
      return -995;
  }

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
    else if (strcmp(argv[argn], "--help") == 0 ||
             strcmp(argv[argn], "-h") == 0)
      return usage(argv[0]);
    else if (strcmp(argv[argn], "--param_card") == 0) {
      param_card = argv[argn + 1];
      argn++;
    }
    else if (strcmp(argv[argn], "--json_file") == 0) {
      json_file = argv[argn + 1];
      json_file_bool = true;
      argn++;
    }
    else if (strcmp(argv[argn], "--device_info") == 0) {
        device_info = true;
    }
    else if (strcmp(argv[argn], "--device_id") == 0) {
        d_id_str = argv[argn + 1];
        if (check_digits(d_id_str)) {
            d_id = std::stoi(d_id_str);
            device_chosen = true;
            if (d_id >= devices.size()) {
                std::cout << "Invalid device_id. Please choose device from: " << std::endl;
                for (int i=0; i < devices.size(); i++) {
                    const auto &device = devices[i];
                    auto d_name = device.get_info<sycl::info::device::name>();
                    std::cout << "    " << i << ": " << d_name << std::endl;
                }
                std::cout << "Terminating. Exit Code: -999" << std::endl;
                return -999;
            }
        } else {
            std::cout << std::endl << "Invalid device_id, must be integer. Terminating. Exit Code: -998" << std::endl << std::endl;
            return usage(argv[0], -998);
        }
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

  if (device_info) {
      if (device_chosen) {
          print_device_info(devices[d_id]);
          std::cout << "Terminating. Exit Code: -997" << std::endl;
          return -997;
      } else {
          for (int i=0; i < devices.size(); i++) {
              const auto &device = devices[i];
              std::cout << "device_id " << i << ":" << std::endl;
              print_device_info(devices[i]);
          }
          std::cout << "Terminating. Exit Code: -996" << std::endl;
          return -996;
      }
  }

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

  sycl::queue q_ct1 = sycl::queue(devices[d_id]);

  auto device = q_ct1.get_device();
  std::cout << "Selected " << device.get_info<sycl::info::device::name>()
            << " on platform "
            << device.get_info<sycl::info::device::platform>().get_info<sycl::info::platform::name>()
            << std::endl;

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
  process.initProc(param_card);
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
  auto hstMomenta   = hstMakeUnique<fptype>( nMomenta ); // AOSOA[npagM][npar][np4][neppM] (previously was: lp)
  auto hstIsGoodHel = hstMakeUnique<bool  >( ncomb );
  auto hstWeights   = hstMakeUnique<fptype>( nWeights ); // (previously was: meHostPtr)
  auto hstMEs       = hstMakeUnique<fptype>( nMEs ); // (previously was: meHostPtr)

  auto devRnarray   = malloc_device<fptype>( nRnarray, q_ct1 ); // AOSOA[npagR][nparf][np4][neppR] (NB: nevt=npagR*neppR)
  auto devMomenta   = malloc_device<fptype>( nMomenta, q_ct1 ); // (previously was: allMomenta)
  auto devIsGoodHel = malloc_device<bool  >( ncomb, q_ct1 );
  auto devWeights   = malloc_device<fptype>( nWeights, q_ct1 ); // (previously was: meDevPtr)
  auto devMEs       = malloc_device<fptype>( nMEs, q_ct1 ); // (previously was: meDevPtr)
  auto devcHel      = malloc_device<int   >( ncomb*npar, q_ct1 );
  auto devcIPC      = malloc_device<fptype>( 6, q_ct1 );
  auto devcIPD      = malloc_device<fptype>( 2, q_ct1 );
  //auto tmp_cHel = process.get_cHel();
  q_ct1.memcpy( devcHel, process.get_cHel_ptr(), ncomb*npar*sizeof(int) ).wait();
  q_ct1.memcpy( devcIPC, process.get_cIPC_ptr(), 6*sizeof(fptype) ).wait();
  q_ct1.memcpy( devcIPD, process.get_cIPD_ptr(), 2*sizeof(fptype) ).wait();
  auto devcNGoodHel = malloc_device<int   >( 1, q_ct1 ); 
  auto devcGoodHel  = malloc_device<int   >( ncomb, q_ct1 ); 

#if defined MGONGPU_CURAND_ONHOST or defined MGONGPU_COMMONRAND_ONHOST
  const int nbytesRnarray = nRnarray * sizeof(fptype);
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
    std::vector<double> commonRnd = commonRandomPromises[iiter].get_future().get();
    assert( nRnarray == static_cast<int>( commonRnd.size() ) );
    // NB (PR #45): memcpy is strictly needed only in CUDA (copy to pinned memory), but keep it also in C++ for consistency
#ifdef SYCL_LANGUAGE_VERSION
    q_ct1.memcpy(devRnarray, commonRnd.data(), nbytesRnarray).wait();
#else
    memcpy( hstRnarray.get(), commonRnd.data(), nRnarray * sizeof(hstRnarray[0]) );
#endif
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
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
              rambo2toNm0::getMomentaInitial(energy, devMomenta, item_ct1);
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
      q_ct1.submit([&](sycl::handler &cgh) {
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
              rambo2toNm0::getMomentaFinal(energy, devRnarray,
                                            devMomenta,
                                            devWeights, item_ct1);
            });
      });
    q_ct1.wait();

    }
#else
    rambo2toNm0::getMomentaFinal( energy, hstRnarray.get(), hstMomenta.get(), hstWeights.get(), nevt );
#endif
    //std::cout << "Got final momenta" << std::endl;

#ifdef SYCL_LANGUAGE_VERSION
    // --- 2c. CopyDToH Weights
    const std::string cwgtKey = "2c CpDTHwgt";
    rambtime += timermap.start( cwgtKey );
    q_ct1.memcpy(hstWeights.get(), devWeights, nbytesWeights).wait();

    // --- 2d. CopyDToH Momenta
    const std::string cmomKey = "2d CpDTHmom";
    rambtime += timermap.start( cmomKey );
    q_ct1.memcpy(hstMomenta.get(), devMomenta, nbytesMomenta).wait();
#endif

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
          cgh.parallel_for(
              sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                    sycl::range<3>(1, 1, gputhreads),
                                sycl::range<3>(1, 1, gputhreads)),
              [=](sycl::nd_item<3> item_ct1) {
                    Proc::sigmaKin_getGoodHel(devMomenta, devMEs, devIsGoodHel, item_ct1, devcHel, devcIPC, devcIPD);
              });
        });
      q_ct1.wait();
      }
      // ... 0d2. Copy back good helicity mask to the host
      q_ct1.memcpy(hstIsGoodHel.get(), devIsGoodHel, nbytesIsGoodHel).wait();
      // ... 0d3. Copy back good helicity list to constant memory on the device
          q_ct1.memset( devcGoodHel, 0, ncomb*sizeof(int) ).wait();
	  {
	      q_ct1.single_task([=] {
                  Proc::sigmaKin_setGoodHel(devIsGoodHel, devcNGoodHel, devcGoodHel);
	      });
              q_ct1.wait();
	  }
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
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                  sycl::range<3>(1, 1, gputhreads),
                              sycl::range<3>(1, 1, gputhreads)),
            [=](sycl::nd_item<3> item_ct1) {
                Proc::sigmaKin(
                devMomenta, devMEs, item_ct1,

		devcHel, devcIPC, devcIPD,
		devcNGoodHel, devcGoodHel);
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
    q_ct1.memcpy(hstMEs.get(), devMEs, nbytesMEs).wait();
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
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 0*neppM + ieppM] // AOSOA[ipagM][ipar][0][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 1*neppM + ieppM] // AOSOA[ipagM][ipar][1][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 2*neppM + ieppM] // AOSOA[ipagM][ipar][2][ieppM]
                    << std::setw(14) << hstMomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + 3*neppM + ieppM] // AOSOA[ipagM][ipar][3][ieppM]
                    << std::endl
                    << std::defaultfloat; // default format: affects all floats
        }
        std::cout << std::string(80, '-') << std::endl;
        // Display matrix elements
        std::cout << " Matrix element = "
                  << hstMEs[ievt] << " GeV^" << meGeVexponent << std::endl; // FIXME: assume process.nprocesses == 1
        std::cout << std::string(80, '-') << std::endl;
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
  double stdwtim = std::sqrt( sqswtim / niter - meanwtim * meanwtim );

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
#ifdef SYCL_LANGUAGE_VERSION
  // --- 9aa. Free device memory 
  const std::string syclfrKey = "9aa sycl_free";
  timermap.start( syclfrKey );

  sycl::free( devRnarray   , q_ct1);
  sycl::free( devMomenta   , q_ct1);
  sycl::free( devIsGoodHel , q_ct1);
  sycl::free( devWeights   , q_ct1);
  sycl::free( devMEs       , q_ct1);
  sycl::free( devcHel      , q_ct1);
  sycl::free( devcIPC      , q_ct1);
  sycl::free( devcIPD      , q_ct1);
#endif
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
#if defined MGONGPU_CURAND_ONDEVICE
             << "\"DEVICE (CUDA code)\"," << std::endl
#elif defined MGONGPU_CURAND_ONHOST
             << "\"HOST (CUDA code)\"," << std::endl
#endif
#else
             << "\"HOST (C++ code)\"," << std::endl
#endif
             << "\"NumberOfEntries\": " << niter << "," << std::endl
      //<< std::scientific // Not sure about this
             << "\"TotalTimeInWaveFuncs\": "
             << "\"" << std::to_string(sumwtim) << " sec\"," << std::endl
             << "\"MeanTimeInWaveFuncs\": "
             << "\"" << std::to_string(meanwtim) << " sec\"," << std::endl
             << "\"StdDevTimeInWaveFuncs\": "
             << "\"" << std::to_string(stdwtim) << " sec\"," << std::endl
             << "\"MinTimeInWaveFuncs\": "
             << "\"" << std::to_string(minwtim) << " sec\"," << std::endl
             << "\"MaxTimeInWaveFuncs\": "
             << "\"" << std::to_string(maxwtim) << " sec\"," << std::endl
      //<< "ProcessID:                = " << getpid() << std::endl
      //<< "NProcesses                = " << process.nprocesses << std::endl
             << "\"TotalEventsComputed\": " << nevtALL << "," << std::endl
             << "\"RamboEventsPerSec\": "
             << "\"" << std::to_string(nevtALL/sumrtim) << " sec^-1\"," << std::endl
             << "\"MatrixElemEventsPerSec\": "
             << "\"" << std::to_string(nevtALL/sumwtim) << " sec^-1\"," << std::endl
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
