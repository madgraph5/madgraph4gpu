//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"

#include "CPPProcess.h"

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ WEIGHTED<=2 @1

#if defined(ALPAKA) || defined(__CUDACC__)
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
  //constexpr fptype cIPC[4] = { ... };
  //constexpr fptype cIPD[2] = { ... };
#ifdef ALPAKA
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPC[4];
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPD[2];
#elif defined(__CUDACC__)
  __device__ __constant__ fptype cIPC[4];
  __device__ __constant__ fptype cIPD[2];
#else
  static fptype cIPC[4];
  static fptype cIPD[2];
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef ALPAKA
  ALPAKA_STATIC_ACC_MEM_CONSTANT short cHel[ncomb][npar];
  ALPAKA_STATIC_ACC_MEM_CONSTANT int cNGoodHel;
  ALPAKA_STATIC_ACC_MEM_CONSTANT int cGoodHel[ncomb];
#elif defined(__CUDACC__)
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ int cNGoodHel; // FIXME: assume process.nprocesses == 1 for the moment (eventually cNGoodHel[nprocesses]?)
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
  static int cNGoodHel; // FIXME: assume process.nprocesses == 1 for the moment (eventually cNGoodHel[nprocesses]?)
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef ALPAKA
  template< typename T_Acc>
  ALPAKA_FN_ACC
  void calculate_wavefunctions( T_Acc const &acc,
                                int ihel,
                                const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#else
  __device__
  INLINE
  void calculate_wavefunctions( int ihel,
                                const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                                , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
#endif
                                )
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: nevt %d\n", nevt );
#endif

    // The number of colors
    constexpr int ncolor = 2;

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1]; // invariant amplitude for one given Feynman diagram

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    cxtype_sv jamp_sv[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===
#if !defined(ALPAKA) && !defined(__CUDACC__)
    const int npagV = nevt / neppV;
    // ** START LOOP ON IPAGV **
#ifdef _OPENMP
    // (NB gcc9 or higher, or clang, is required)
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default(none) shared(allmomenta,allMEs,cHel,cIPC,cIPD,ihel,npagV) private (amp_sv,w_sv,jamp_sv)
#endif
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif
    {
      // Reset color flows (reset jamp_sv) at the beginning of a new event or event page
      for( int i=0; i<ncolor; i++ ){ jamp_sv[i] = cxzero_sv(); }

      // *** DIAGRAM 1 OF 3 ***

      // Wavefunction(s) for diagram number 1
#ifdef ALPAKA
      vxxxxx(acc, allmomenta, 0., cHel[ihel][0], -1, w_sv[0], 0);
#elif defined(__CUDACC__)
      vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w_sv[0], 0);
#else
      vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w_sv[0], ipagV, 0);
#endif

#ifdef ALPAKA
      vxxxxx(acc, allmomenta, 0., cHel[ihel][1], -1, w_sv[1], 1);
#elif defined(__CUDACC__)
      vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w_sv[1], 1);
#else
      vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w_sv[1], ipagV, 1);
#endif

#ifdef ALPAKA
      oxxxxx(acc, allmomenta, cIPD[0], cHel[ihel][2], +1, w_sv[2], 2);
#elif defined(__CUDACC__)
      oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w_sv[2], 2);
#else
      oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w_sv[2], ipagV, 2);
#endif

#ifdef ALPAKA
      ixxxxx(acc, allmomenta, cIPD[0], cHel[ihel][3], -1, w_sv[3], 3);
#elif defined(__CUDACC__)
      ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w_sv[3], 3);
#else
      ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w_sv[3], ipagV, 3);
#endif

#ifdef ALPAKA
      VVV1P0_1( acc, w_sv[0], w_sv[1], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_sv[4] );
#else
      VVV1P0_1( w_sv[0], w_sv[1], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_sv[4] );
#endif

      // Amplitude(s) for diagram number 1
#ifdef ALPAKA
      FFV1_0( acc, w_sv[3], w_sv[2], w_sv[4], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#else
      FFV1_0( w_sv[3], w_sv[2], w_sv[4], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#endif
      jamp_sv[0] += +cxtype(0,1)*amp_sv[0];
      jamp_sv[1] -= cxtype(0,1)*amp_sv[0];

      // *** DIAGRAM 2 OF 3 ***

      // Wavefunction(s) for diagram number 2
#ifdef ALPAKA
      FFV1_1( acc, w_sv[2], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cIPD[0], cIPD[1], w_sv[4] );
#else
      FFV1_1( w_sv[2], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cIPD[0], cIPD[1], w_sv[4] );
#endif

      // Amplitude(s) for diagram number 2
#ifdef ALPAKA
      FFV1_0( acc, w_sv[3], w_sv[4], w_sv[1], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#else
      FFV1_0( w_sv[3], w_sv[4], w_sv[1], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 3 OF 3 ***

      // Wavefunction(s) for diagram number 3
#ifdef ALPAKA
      FFV1_2( acc, w_sv[3], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cIPD[0], cIPD[1], w_sv[4] );
#else
      FFV1_2( w_sv[3], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cIPD[0], cIPD[1], w_sv[4] );
#endif

      // Amplitude(s) for diagram number 3
#ifdef ALPAKA
      FFV1_0( acc, w_sv[4], w_sv[2], w_sv[1], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#else
      FFV1_0( w_sv[4], w_sv[2], w_sv[1], cxmake( cIPC[2], cIPC[3] ), &amp_sv[0] );
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called CPPProcess::matrix_1_gg_ttx()?)

      // The color matrix
      static constexpr fptype denom[ncolor] = {3, 3};
      static constexpr fptype cf[ncolor][ncolor] = {
      {16, -2},
      {-2, 16}};

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros
      for( int icol = 0; icol < ncolor; icol++ )
      {
        cxtype_sv ztemp_sv = cxzero_sv();
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp_sv += cf[icol][jcol] * jamp_sv[jcol];
        deltaMEs += cxreal( ztemp_sv * cxconj( jamp_sv[icol] ) ) / denom[icol];
      }

      // *** STORE THE RESULTS ***

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#if defined(ALPAKA) || defined(__CUDACC__)
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      allMEs[ievt] += deltaMEs;
      //printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );
#else
      allMEs[ipagV] += deltaMEs;
      //printf( "calculate_wavefunction: %6d %2d %f\n", ipagV, ihel, allMEs[ipagV] ); // FIXME for MGONGPU_CPPSIMD
#endif
    }

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( int numiterations,
                          int ngpublocks,
                          int ngputhreads,
                          bool verbose,
                          bool debug )
    : m_numiterations( numiterations )
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    , m_pars( 0 )
    , m_masses()
  {
    // Helicities for the process - nodim
    static constexpr short tHel[ncomb][mgOnGpu::npar] = {
      {-1, -1, -1, -1},
      {-1, -1, -1, 1},
      {-1, -1, 1, -1},
      {-1, -1, 1, 1},
      {-1, 1, -1, -1},
      {-1, 1, -1, 1},
      {-1, 1, 1, -1},
      {-1, 1, 1, 1},
      {1, -1, -1, -1},
      {1, -1, -1, 1},
      {1, -1, 1, -1},
      {1, -1, 1, 1},
      {1, 1, -1, -1},
      {1, 1, -1, 1},
      {1, 1, 1, -1},
      {1, 1, 1, 1}};
#ifdef ALPAKA
    auto extent = alpaka::Vec<alpaka::DimInt<1u>, int>{ncomb * mgOnGpu::npar};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcHel = alpaka::createStaticDevMemView((short*)cHel, device, extent);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, short, alpaka::DimInt<1u>, int> bufHostHel((short*)tHel, host, extent);
    alpaka::memcpy(streamObject, viewcHel, bufHostHel, extent);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * mgOnGpu::npar * sizeof(short) ) );
#else
    memcpy( cHel, tHel, ncomb * mgOnGpu::npar * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#if !defined(ALPAKA) && !defined(__CUDACC__)
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    assert( neppV == mgOnGpu::neppM );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  // Initialize process
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
    }
    m_pars->setDependentParameters();
    m_pars->setDependentCouplings();

    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->mdl_MT );
    m_masses.push_back( m_pars->mdl_MT );

    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const cxtype tIPC[2] = { cxmake( m_pars->GC_10 ), cxmake( m_pars->GC_11 ) };
    const fptype tIPD[2] = { (fptype)m_pars->mdl_MT, (fptype)m_pars->mdl_WT };
#ifdef ALPAKA
    auto extentcx3 = alpaka::Vec<alpaka::DimInt<1u>, int>{6};
    auto extentfp2 = alpaka::Vec<alpaka::DimInt<1u>, int>{2};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcIPC = alpaka::createStaticDevMemView((fptype*)cIPC, device, extentcx3);
    auto viewcIPD = alpaka::createStaticDevMemView((fptype*)cIPD, device, extentfp2);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, fptype, alpaka::DimInt<1u>, int> bufHostIPC((fptype *)tIPC, host, extentcx3);
    alpaka::ViewPlainPtr<cupla::AccHost, fptype, alpaka::DimInt<1u>, int> bufHostIPD((fptype *)tIPD, host, extentfp2);
    alpaka::memcpy(streamObject, viewcIPC, bufHostIPC, extentcx3);
    alpaka::wait( streamObject );
    alpaka::memcpy(streamObject, viewcIPD, bufHostIPD, extentfp2);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 2 * sizeof(cxtype) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(fptype) ) );
#else
    memcpy( cIPC, tIPC, 2 * sizeof(cxtype) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
#endif

    //std::cout << std::setprecision(17) << "tIPC[0] = " << tIPC[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[1] = " << tIPC[1] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[2] = " << tIPC[2] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[0] = " << tIPD[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[1] = " << tIPD[1] << std::endl;
  }

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler()
  {
    std::stringstream out;
    // CUDA version (NVCC)
#ifdef __CUDACC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
    out << " (";
#endif
    // ICX version (either as CXX or as host compiler inside NVCC)
#if defined __INTEL_COMPILER
#error "icc is no longer supported: please use icx"
#elif defined __INTEL_LLVM_COMPILER // alternative: __INTEL_CLANG_COMPILER
    out << "icx " << __INTEL_LLVM_COMPILER << " (";
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype(&pclose)> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if ( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while ( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
#endif
#else
    out << "clang UNKNOWKN";
#endif
#else
    // GCC version (either as CXX or as host compiler inside NVCC)
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
#ifdef ALPAKA
    out << " (ALPAKA)";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

#if defined(ALPAKA) || defined(__CUDACC__)
#ifdef ALPAKA
template< typename T_Acc >
ALPAKA_FN_ACC
  void sigmaKin_getGoodHel::operator()( T_Acc const &acc,
                            const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel ) const      // output: isGoodHel[ncomb] - device array
#else
  __global__
  void sigmaKin_getGoodHel( 
                            const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel )            // output: isGoodHel[ncomb] - device array
#endif
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef ALPAKA
      calculate_wavefunctions( acc, ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#endif
      if ( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                            , const int nevt )           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int maxtry0 = ( neppV > 16 ? neppV : 16 ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    fptype_sv allMEsLast[maxtry0/neppV] = { 0 };
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for ( int ipagV = 0; ipagV < maxtry/neppV; ++ipagV )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ipagV] = fptype_sv{0}; // all zeros
    }
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allMEs, maxtry );
      for ( int ipagV = 0; ipagV < maxtry/neppV; ++ipagV )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        const bool differs = maskor( allMEs[ipagV] != allMEsLast[ipagV] ); // true if any of the neppV events differs
        if ( differs )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ipagV] = allMEs[ipagV]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel = 0; // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        //goodHel[nGoodHel[0]] = ihel; // FIXME: assume process.nprocesses == 1 for the moment
        //nGoodHel[0]++; // FIXME: assume process.nprocesses == 1 for the moment
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef ALPAKA
    auto extent1 = alpaka::Vec<alpaka::DimInt<1u>, int>{1};
    auto extentn = alpaka::Vec<alpaka::DimInt<1u>, int>{ncomb};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcNgood = alpaka::createStaticDevMemView(&cNGoodHel, device, extent1);
    auto viewcGood = alpaka::createStaticDevMemView((int*)cGoodHel, device, extentn);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, int, alpaka::DimInt<1u>, int> bufHostNgood(&nGoodHel, host, extent1);
    alpaka::ViewPlainPtr<cupla::AccHost, int, alpaka::DimInt<1u>, int> bufHostGood((int*)goodHel, host, extentn);
    alpaka::memcpy(streamObject, viewcNgood, bufHostNgood, extent1);
    alpaka::wait( streamObject );
    alpaka::memcpy(streamObject, viewcGood, bufHostGood, extentn);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) ); // FIXME: assume process.nprocesses == 1 for the moment
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb*sizeof(int) ) );
#else
    cNGoodHel = nGoodHel;
    for ( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC
  void sigmaKin::operator()( T_Acc const &acc,
                             const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                             fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                           ) const
#else
  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
#endif
  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    const int denominators = 256; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

#if defined(ALPAKA) || defined(__CUDACC__)
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "sigmakin: ievt %d\n", ievt );
#endif

    // Start sigmaKin_lines
    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#if defined (ALPAKA) || defined( __CUDACC__)
    allMEs[ievt] = 0;
#else
    const int npagV = nevt/neppV;
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      allMEs[ipagV] = fptype_sv{ 0 };
    }
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef ALPAKA
      calculate_wavefunctions( acc, ihel, allmomenta, allMEs );
#elif defined(__CUDACC__)
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs, nevt );
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#if defined(ALPAKA) || defined(__CUDACC__)
    allMEs[ievt] /= (fptype)denominators;
#else
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      allMEs[ipagV] /= (fptype)denominators;
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================

// This was initially added to both C++ and CUDA in order to avoid RDC in CUDA (issue #51)
// This is now also needed by C++ LTO-like optimizations via inlining (issue #229)
#include "HelAmps_sm.cc"

//==========================================================================

#ifdef ALPAKA
template void gProc::sigmaKin_getGoodHel::operator()<cupla::Acc>(cupla::Acc const&, const fptype_sv* , fptype_sv*, bool *) const;
template void gProc::sigmaKin::operator()<cupla::Acc>(cupla::Acc const&, const fptype_sv* , fptype_sv*) const;
#endif

