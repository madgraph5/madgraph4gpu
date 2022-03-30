#ifndef MGONGPUBRIDGE_H
#define MGONGPUBRIDGE_H 1
// Includes from mg4gpu code 
#include <CL/sycl.hpp>
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

#include "read_slha.h"
#include "Parameters_sm.h"
#include "CPPProcess.h"

using mgOnGpu::np4;
using mgOnGpu::nparf;
using mgOnGpu::npar;
using mgOnGpu::ncomb; // Number of helicity combinations

// those should become fortran parameters passed in here
const int gpublocks = 1;   // 1024;
const int gputhreads = 16; // 256; NB_PAGE???

// helper functions

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

/**
const int evnt_n = 4;  // the number of events
const int part_n = 4;  // number of in/out particles inside an event
const int mome_n = 3;  // number of momenta of one particle (usually 4)
const int strd_n = 2;  // stride length for aosoa data (# adjacent events)
const int array_bytes = evnt_n * part_n * mome_n * sizeof(T);
*/
template <typename T>
void dev_transpose(const T *in, T *out, size_t pos, const int evt) {
  constexpr int part = mgOnGpu::npar;
  constexpr int mome = mgOnGpu::np4;
  constexpr int strd = mgOnGpu::neppM;
  int arrlen = evt * part * mome;

  if (pos < arrlen) {
    int page_i = pos / (strd * mome * part);
    int rest_1 = pos % (strd * mome * part);
    int part_i = rest_1 / (strd * mome);
    int rest_2 = rest_1 % (strd * mome);
    int mome_i = rest_2 / strd;
    int strd_i = rest_2 % strd;
    int inpos = (page_i * strd + strd_i) // event number
                    * (part * mome)      // event size (pos of event)
                + part_i * mome          // particle inside event
                + mome_i;                // momentum inside particle

    out[pos] = in[inpos];
  }
}

// *****************************************************************************

/**
 * A templated class for calling the SYCL C++ matrix element calculations of
 * the event generation workflow. The template parameter is used for the
 * precision of the calculations (float or double)
 *
 * The fortran momenta passed in are in the form of
 *   DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
 * where the dimensions are <# momenta>, <# of particles>, <# events>
 */
template<class T>
class Bridge {
public:
    /**
     * class constructor
     *
     * @param evt number of events (NB_PAGE, vector.inc)
     * @param par number of particles / event (NEXTERNAL, nexternal.inc)
     * @param mom number of momenta / particle
     * @param str stride length
     * @param ncomb number of good helicities (ncomb, mgOnGpuConfig.h)
     */
    Bridge();
    ~Bridge();
  
    /**
     * sequence to be executed for matrix element calculation
     *
     * @param momenta memory address of the input 4-momenta
     * @param mes memory address of the output matrix elements
     */
    void generate_matrix_elements(T *momenta, T *mes);


private:
    bool m_goodHelsCalculated = false; ///< have the good helicities been calculated?

    uint32_t d_id = 0; // default to device 0 FIXME need to set device

    bool device_chosen = true;
    bool device_info = false;
    std::vector<sycl::device> devices = sycl::device::get_devices();

    //std::string param_card = "../Cards/param_card.dat"; //FIXME need to set this location
    std::string param_card = "/home/nnichols/mg_nsn_fork/epochX/sycl/gg_ttgg.auto/Cards/param_card.dat";
    
    const int neppR = mgOnGpu::neppR; // ASA layout: constant at compile-time
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time

    const int ndim = gpublocks * gputhreads; // number of threads in one GPU grid
    const int nevt = ndim; // number of events in one iteration == number of GPU threads

    sycl::queue q_ct1;

    // set memory
    const unsigned int nMomenta = np4*npar*nevt; // (NB: nevt=npagM*neppM for ASA layouts)
    const unsigned int nMEs     = nevt;

    fptype *devMEs;
    short *devcHel;
    fptype *devcIPC;
    fptype *devcIPD;
    int *devcGoodHel;
    int *devcNGoodHel;
    fptype *devMomenta;
    bool *devIsGoodHel;

    const int nbytesMomenta = nMomenta * sizeof(fptype);
    const int nbytesIsGoodHel = ncomb * sizeof(bool);
    const int nbytesMEs = nMEs * sizeof(fptype);
};

// *****************************************************************************
// Implementations of class Bridge member functions
//
// *****************************************************************************

template<class T>
Bridge<T>::Bridge() {

    //Print device info
    for (int i=0; i < devices.size(); i++) {
        const auto &device = devices[i];
        std::cout << "device_id " << i << ":" << std::endl;
        print_device_info(devices[i]);
    }

    sycl::queue q_ct1 = sycl::queue(devices[d_id]);

    devMEs          = sycl::malloc_device<fptype>( nMEs, q_ct1 ); // (previously was: meDevPtr)
    devcHel         = sycl::malloc_device<short   >( ncomb*npar, q_ct1 );
    devcIPC         = sycl::malloc_device<fptype>( 6, q_ct1 );
    devcIPD         = sycl::malloc_device<fptype>( 2, q_ct1 );
    devMomenta      = sycl::malloc_device<fptype>( nMomenta, q_ct1 ); // (previously was: allMomenta)
    devMomentaF     = sycl::malloc_device<fptype>( nMomenta, q_ct1 ); // (previously was: allMomenta)
    devcGoodHel     = sycl::malloc_device<int   >( ncomb, q_ct1 ); 
    devIsGoodHel    = sycl::malloc_device<bool  >( ncomb, q_ct1 );
    devcNGoodHel    = sycl::malloc_device<int   >( 1, q_ct1 ); 

    // Create a process object
    Proc::CPPProcess process( 1, gpublocks, gputhreads, false );

    // Read param_card and set parameters
    process.initProc(param_card);

    //Copy data to device
    q_ct1.memcpy( devcHel, process.get_cHel_ptr(), ncomb*npar*sizeof(short) ).wait();
    q_ct1.memcpy( devcIPC, process.get_cIPC_ptr(), 6*sizeof(fptype) ).wait();
    q_ct1.memcpy( devcIPD, process.get_cIPD_ptr(), 2*sizeof(fptype) ).wait();
}

template<class T>
void Bridge<T>::generate_matrix_elements(T *momenta, T *mes) { 
    auto _devMomenta = devMomenta;
    auto _devMomentaF = devMomentaF;
    auto _devMEs = devMEs;
    auto _devIsGoodHel = devIsGoodHel;
    auto _devcNGoodHel = devcNGoodHel;
    auto _devcGoodHel = devcGoodHel;
    auto _devcHel = devcHel;
    auto _devcIPC = devcIPC;
    auto _devcIPD = devcIPD;
    // Copy momenta to device
    q_ct1.memcpy(_devMomentaF, momenta, nbytesMomenta).wait(); //FIXME might need to transpose

    cgh.parallel_for(sycl::range<1>(nMomenta),
            [=](sycl::id<1> item) {
        size_t ievt = item.get(0);
        dev_transpose(_devMomentaF, _devMomenta, ievt, nevt);
    });

    // Set good helicities
    if (!m_goodHelsCalculated) {
        // ... 0a1. Compute good helicity mask on the device
        q_ct1.submit([&](sycl::handler &cgh) {
            cgh.parallel_for(sycl::range<1>(ndim),
                    [=](sycl::id<1> item) {
                size_t ievt = item.get(0);
                Proc::sigmaKin_getGoodHel(_devMomenta, _devMEs, _devIsGoodHel, ievt, _devcHel, _devcIPC, _devcIPD);
            });
        });
        q_ct1.wait();
        // ... 0a2. Copy back good helicity list to constant memory on the device
        q_ct1.memset( devcGoodHel, 0, ncomb*sizeof(int) ).wait();
        q_ct1.submit([&](sycl::handler &cgh) {
            cgh.single_task(
                    [=] (sycl::kernel_handler cgh) {
                Proc::sigmaKin_setGoodHel(_devIsGoodHel, _devcNGoodHel, _devcGoodHel);
            });
        });
        q_ct1.wait();
        m_goodHelsCalculated = true;
    }

    // --- 2a. generate MEs (sigmaKin)
    q_ct1.submit([&](sycl::handler &cgh) {
        cgh.parallel_for(sycl::range<1>(ndim),
                [=](sycl::id<1> item) {
            size_t ievt = item.get(0);
            Proc::sigmaKin( _devMomenta, _devMEs, ievt, _devcHel, _devcIPC, _devcIPD,
                    _devcNGoodHel, _devcGoodHel);
        });
    });
    q_ct1.wait();

    // --- 2b. CopyDToH MEs
    q_ct1.memcpy(mes, _devMEs, nbytesMEs).wait();

}

template<class T>
Bridge<T>::~Bridge() { 
    sycl::free( devMomenta   , q_ct1);
    sycl::free( devIsGoodHel , q_ct1);
    sycl::free( devMEs       , q_ct1);
    sycl::free( devcHel      , q_ct1);
    sycl::free( devcIPC      , q_ct1);
    sycl::free( devcIPD      , q_ct1);
    sycl::free( devcNGoodHel , q_ct1);
    sycl::free( devcGoodHel  , q_ct1);
}

template<class T>
Bridge<T> make_bridge() { return Bridge<T>();}
#endif // MGONGPUBRIDGE_H
