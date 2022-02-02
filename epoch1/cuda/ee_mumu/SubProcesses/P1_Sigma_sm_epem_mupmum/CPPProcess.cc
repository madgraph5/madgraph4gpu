//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"
#include "MemoryAccess.h"

#include "CPPProcess.h"

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4; // 4: the dimension of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // number of particles in total (initial + final)
  using mgOnGpu::ncomb; // number of helicity combinations

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
#ifdef MGONGPU_HARDCODE_CIPC
  __device__ const fptype cIPC[6] = { 0, -0.30795376724436879, 0, -0.28804415396362731, 0, 0.082309883272248419 };
  __device__ const fptype cIPD[2] = { 91.188000000000002, 2.4414039999999999 };
#else
#ifdef __CUDACC__
  __device__ __constant__ fptype cIPC[6];
  __device__ __constant__ fptype cIPD[2];
#else
  static fptype cIPC[6];
  static fptype cIPD[2];
#endif
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef __CUDACC__
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
  // NB: calculate_wavefunctions ADDS |M|^2 for given ihel to running sum of |M|^2 over helicities for given event(s)
  __device__
  INLINE
  void calculate_wavefunctions( int ihel,
                                const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM], nevt=npagM*neppM
                                fptype* allMEs            // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifndef __CUDACC__
                                , const int nevt          // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
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
    const int ncolor = 1;

    // The color matrix
    const fptype denom[ncolor] = {1};
    const fptype cf[ncolor][ncolor] = {{1}};

    // Local variables for the given CUDA event (ievt)
    // Local variables for the given C++ event page (ipagV)
    cxtype_sv w_sv[nwf][nw6]; // w_v[5][6]
    cxtype_sv amp_sv[1]; // was 2

#ifndef __CUDACC__
    const int npagV = nevt / neppV;
#ifdef MGONGPU_CPPSIMD
    const bool isAligned_allMEs = ( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // require SIMD-friendly alignment by at least neppV*sizeof(fptype)
#endif
    // ** START LOOP ON IPAGV **
#ifdef _OPENMP
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#if not defined __clang__ && defined __GNUC__ && __GNUC__ < 9
#pragma omp parallel for default(none) \
  shared(allmomenta,allMEs,cf,cHel,cIPC,cIPD,denom,ihel,cNGoodHel) private (amp_sv,w_sv)
#elif defined MGONGPU_CPPSIMD
#pragma omp parallel for default(none) \
  shared(allmomenta,allMEs,cf,cHel,cIPC,cIPD,denom,ihel,cNGoodHel,npagV,isAligned_allMEs) private (amp_sv,w_sv)
#else
#pragma omp parallel for default(none) \
  shared(allmomenta,allMEs,cf,cHel,cIPC,cIPD,denom,ihel,cNGoodHel,npagV) private (amp_sv,w_sv)
#endif
#endif
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif
    {
#ifdef __CUDACC__
#ifndef MGONGPU_TEST_DIVERGENCE
      // NB: opzxxx only reads pz (not E,px,py)
      opzxxx( allmomenta, cHel[ihel][0], -1, w_sv[0], 0 );
      //oxxxxx( allmomenta, 0, cHel[ihel][0], -1, w_sv[0], 0 ); // tested ok (much slower)
#else
      if ( ( blockDim.x * blockIdx.x + threadIdx.x ) % 2 == 0 )
        opzxxx( allmomenta, cHel[ihel][0], -1, w_sv[0], 0 );
      else
        oxxxxx( allmomenta, 0, cHel[ihel][0], -1, w_sv[0], 0 ); // tested ok (much slower)
#endif
#else
      opzxxx( p4IparIpagV( allmomenta, 0, ipagV ), cHel[ihel][0], -1, w_sv[0] );
      //oxxxxx( p4IparIpagV( allmomenta, 0, ipagV ), 0, cHel[ihel][0], -1, w_sv[0] ); // tested ok (slower)
#endif

#ifdef __CUDACC__
      // NB: imzxxx only reads pz (not E,px,py)
      imzxxx( allmomenta, cHel[ihel][1], +1, w_sv[1], 1 );
      //ixxxxx( allmomenta, 0, cHel[ihel][1], +1, w_sv[1], 1 ); // tested ok (slower)
#else
      imzxxx( p4IparIpagV( allmomenta, 1, ipagV ), cHel[ihel][1], +1, w_sv[1] );
      //ixxxxx( p4IparIpagV( allmomenta, 1, ipagV ), 0, cHel[ihel][1], +1, w_sv[1] ); // tested ok (a bit slower)
#endif

#ifdef __CUDACC__
      // NB: ixzxxx reads all E,px,py,pz
      ixzxxx( allmomenta, cHel[ihel][2], -1, w_sv[2], 2 );
      //ixxxxx( allmomenta, 0, cHel[ihel][2], -1, w_sv[2], 2 ); // tested ok (a bit slower)
#else
      ixzxxx( p4IparIpagV( allmomenta, 2, ipagV ), cHel[ihel][2], -1, w_sv[2] );
      //ixxxxx( p4IparIpagV( allmomenta, 2, ipagV ), 0, cHel[ihel][2], -1, w_sv[2] ); // tested ok (a bit slower)
#endif

#ifdef __CUDACC__
      // NB: oxzxxx reads all E,px,py,pz
      oxzxxx( allmomenta, cHel[ihel][3], +1, w_sv[3], 3 );
      //oxxxxx( allmomenta, 0, cHel[ihel][3], +1, w_sv[3], 3 ); // tested ok (a bit slower)
#else
      oxzxxx( p4IparIpagV( allmomenta, 3, ipagV ), cHel[ihel][3], +1, w_sv[3] );
      //oxxxxx( p4IparIpagV( allmomenta, 3, ipagV ), 0, cHel[ihel][3], +1, w_sv[3] ); // tested ok (a bit slower)
#endif

      // Local variables for the given CUDA event (ievt)
      // Local variables for the given C++ event page (ipagV)
      cxtype_sv jamp_sv[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams

      // --- START Compute amplitudes for all diagrams ---
      FFV1P0_3( w_sv[1], w_sv[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_sv[4] );
      // Amplitude(s) for diagram number 1
      FFV1_0( w_sv[2], w_sv[3], w_sv[4], cxmake( cIPC[0], cIPC[1] ), &amp_sv[0] );
      jamp_sv[0] -= amp_sv[0];

      FFV2_4_3( w_sv[1], w_sv[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w_sv[4] );
      // Amplitude(s) for diagram number 2
      FFV2_4_0( w_sv[2], w_sv[3], w_sv[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp_sv[0] );
      jamp_sv[0] -= amp_sv[0];
      // --- END   Compute amplitudes for all diagrams ---

      // --- START Color matrix algebra ---      
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
      // --- END   Color matrix algebra ---      

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for given ihel to running sum of |M|^2 over helicities for given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      allMEs[ievt] += deltaMEs;
      //if ( cNGoodHel > 0 ) printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );
#else
#ifdef MGONGPU_CPPSIMD
      if ( isAligned_allMEs )
      {
        *reinterpret_cast<fptype_sv*>( &( allMEs[ipagV*neppV] ) ) += deltaMEs;
      }
      else
      {
        for ( int ieppV=0; ieppV<neppV; ieppV++ )
          allMEs[ipagV*neppV + ieppV] += deltaMEs[ieppV];
      }
      //if ( cNGoodHel > 0 )
      //  for ( int ieppV=0; ieppV<neppV; ieppV++ )
      //    printf( "calculate_wavefunction: %6d %2d %f\n", ipagV*neppV+ieppV, ihel, allMEs[ipagV][ieppV] );
#else
      allMEs[ipagV] += deltaMEs;
      //if ( cNGoodHel > 0 ) printf( "calculate_wavefunction: %6d %2d %f\n", ipagV, ihel, allMEs[ipagV] );
#endif
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
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
    static constexpr short tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(short) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    static_assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#ifndef __CUDACC__
    // SANITY CHECK: check that neppR, neppM and neppV are powers of two (https://stackoverflow.com/a/108360)
    auto ispoweroftwo = []( int n ) { return ( n > 0 ) && !( n & ( n - 1 ) ); };
    static_assert( ispoweroftwo( mgOnGpu::neppR ) );
    static_assert( ispoweroftwo( mgOnGpu::neppM ) );
    static_assert( ispoweroftwo( neppV ) );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

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
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );

#ifndef MGONGPU_HARDCODE_CIPC
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const cxtype tIPC[3] = { cxmake( m_pars->GC_3 ), cxmake( m_pars->GC_50 ), cxmake( m_pars->GC_59 ) };
    const fptype tIPD[2] = { (fptype)m_pars->mdl_MZ, (fptype)m_pars->mdl_WZ };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(cxtype) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(fptype) ) );
#else
    memcpy( cIPC, tIPC, 3 * sizeof(cxtype) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
#endif
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
    return out.str();
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM], nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool* isGoodHel )         // output: isGoodHel[ncomb] - device array
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for given ihel to running sum of |M|^2 over helicities for given event(s)
      calculate_wavefunctions( ihel, allmomenta, allMEs );
      if ( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM], nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool* isGoodHel           // output: isGoodHel[ncomb] - device array
                            , const int nevt )        // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    const int maxtry0 = ( neppV > 16 ? neppV : 16 ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    fptype allMEsLast[maxtry0] = { 0 };
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for ( int ievt = 0; ievt < maxtry; ++ievt )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allMEs, maxtry );
      for ( int ievt = 0; ievt < maxtry; ++ievt )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        const bool differs = ( allMEs[ievt] != allMEsLast[ievt] );
        if ( differs )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ievt] = allMEs[ievt]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    // FIXME: assume process.nprocesses == 1 for the moment
    //int nGoodHel[1] = { 0 };
    int nGoodHel = 0;
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) );
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb*sizeof(int) ) );
#else
    cNGoodHel = nGoodHel;
    for ( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  __global__
  void sigmaKin( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs            // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifndef __CUDACC__
                 , const int nevt          // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
    // start sigmakin_lines

    // Denominators: spins, colors and identical particles
    //const int nprocesses = 1;
    //const int denominators[nprocesses] = { 4 };
    const int denominators = 4;

#ifdef __CUDACC__
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "sigmakin: ievt %d\n", ievt );
#else
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
#endif

    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    {
      allMEs[ievt] = 0; // all zeros
    }
#else
    const int npagV = nevt/neppV;    
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for ( int ieppV=0; ieppV<neppV; ieppV++ )
        allMEs[ipagV*neppV + ieppV] = 0; // all zeros
    }
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef __CUDACC__
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs, nevt );
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    {
      allMEs[ievt] /= (fptype)denominators;
    }
#else
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for ( int ieppV=0; ieppV<neppV; ieppV++ )
        allMEs[ipagV*neppV + ieppV] /= (fptype)denominators;
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================

// This was initially added to both C++ and CUDA in order to avoid RDC in CUDA (issue #51)
// This is now also needed by C++ LTO-like optimizations via inlining (issue #229)
#include "../../src/HelAmps_sm.cc"

//==========================================================================

