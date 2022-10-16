//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"

#include "mgOnGpuConfig.h"

#include "CudaRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#endif

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>

// Test ncu metrics for CUDA thread divergence
#undef MGONGPU_TEST_DIVERGENCE
//#define MGONGPU_TEST_DIVERGENCE 1

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ g g WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  using mgOnGpu::np4;   // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar;  // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QCD)
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on running alphas QCD)

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
#ifdef MGONGPU_HARDCODE_PARAM
  __device__ const fptype cIPD[2] = { (fptype)Parameters_sm::mdl_MT, (fptype)Parameters_sm::mdl_WT };
  __device__ const fptype* cIPC = nullptr; // unused as nicoup=0
#else
#ifdef __CUDACC__
  __device__ __constant__ fptype cIPD[2];
  __device__ __constant__ fptype* cIPC = nullptr; // unused as nicoup=0
#else
  static fptype cIPD[2];
  static fptype* cIPC = nullptr; // unused as nicoup=0
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
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  // (similarly, it also ADDS the numerator and denominator for a given ihel to their running sums over helicities)
  __device__ INLINE void /* clang-format off */
  calculate_wavefunctions( int ihel,
                           const fptype* allmomenta,      // input: momenta[nevt*npar*4]
                           const fptype* allcouplings,    // input: couplings[nevt*ndcoup*2]
                           fptype* allMEs                 // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                           , fptype* allNumerators        // output: multichannel numerators[nevt], running_sum_over_helicities
                           , fptype* allDenominators      // output: multichannel denominators[nevt], running_sum_over_helicities
                           , const unsigned int channelId // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
#ifndef __CUDACC__
                           , const int nevt               // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                           )
  //ALWAYS_INLINE // attributes are not permitted in a function definition
  {
#ifdef __CUDACC__
    using namespace mg5amcGpu;
    using M_ACCESS = DeviceAccessMomenta;         // non-trivial access: buffer includes all events
    using E_ACCESS = DeviceAccessMatrixElements;  // non-trivial access: buffer includes all events
    using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = DeviceAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#else
    using namespace mg5amcCpu;
    using M_ACCESS = HostAccessMomenta;         // non-trivial access: buffer includes all events
    using E_ACCESS = HostAccessMatrixElements;  // non-trivial access: buffer includes all events
    using W_ACCESS = HostAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = HostAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = HostAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#endif /* clang-format on */
    mgDebug( 0, __FUNCTION__ );
    //printf( "calculate_wavefunctions: ihel=%2d\n", ihel );
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: nevt %d\n", nevt );
#endif

    // The number of colors
    constexpr int ncolor = 24;

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
    // ** NB: in other words, amplitudes and wavefunctions still have TRIVIAL ACCESS: there is currently no need
    // ** NB: to have large memory structurs for wavefunctions/amplitudes in all events (no kernel splitting yet)!
    //MemoryBufferWavefunctions w_buffer[nwf]{ neppV };
    cxtype_sv w_sv[nwf][nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
    cxtype_sv amp_sv[1];      // invariant amplitude for one given Feynman diagram

    // Proof of concept for using fptype* in the interface
    fptype* w_fp[nwf];
    for( int iwf = 0; iwf < nwf; iwf++ ) w_fp[iwf] = reinterpret_cast<fptype*>( w_sv[iwf] );
    fptype* amp_fp;
    amp_fp = reinterpret_cast<fptype*>( amp_sv );

    // Local variables for the given CUDA event (ievt) or C++ event page (ipagV)
    // [jamp: sum (for one event or event page) of the invariant amplitudes for all Feynman diagrams in a given color combination]
    cxtype_sv jamp_sv[ncolor] = {}; // all zeros (NB: vector cxtype_v IS initialized to 0, but scalar cxype is NOT, if "= {}" is missing!)

    // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===
#ifndef __CUDACC__
    const int npagV = nevt / neppV;
    // ** START LOOP ON IPAGV **
#ifdef _OPENMP
    // (NB gcc9 or higher, or clang, is required)
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default( none ) shared( allmomenta, allMEs, cHel, allcouplings, cIPC, cIPD, ihel, npagV, amp_fp, w_fp ) private( amp_sv, w_sv, jamp_sv )
#endif // _OPENMP
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif // !__CUDACC__
    {
      constexpr size_t nxcoup = ndcoup + nicoup; // both dependent and independent couplings
      const fptype* allCOUPs[nxcoup];
#ifdef __CUDACC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 186 // e.g. <<warning #186-D: pointless comparison of unsigned integer with zero>>
#endif
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        allCOUPs[idcoup] = CD_ACCESS::idcoupAccessBufferConst( allcouplings, idcoup ); // dependent couplings, vary event-by-event
      for( size_t iicoup = 0; iicoup < nicoup; iicoup++ )
        allCOUPs[ndcoup + iicoup] = CI_ACCESS::iicoupAccessBufferConst( cIPC, iicoup ); // independent couplings, fixed for all events
#ifdef __CUDACC__
#pragma nv_diagnostic pop
      // CUDA kernels take input/output buffers with momenta/MEs for all events
      const fptype* momenta = allmomenta;
      const fptype* COUPs[nxcoup];
      for( size_t ixcoup = 0; ixcoup < nxcoup; ixcoup++ ) COUPs[ixcoup] = allCOUPs[ixcoup];
      fptype* MEs = allMEs;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = allNumerators;
      fptype* denominators = allDenominators;
#endif
#else
      // C++ kernels take input/output buffers with momenta/MEs for one specific event (the first in the current event page)
      const int ievt0 = ipagV * neppV;
      const fptype* momenta = M_ACCESS::ieventAccessRecordConst( allmomenta, ievt0 );
      const fptype* COUPs[nxcoup];
      for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
        COUPs[idcoup] = CD_ACCESS::ieventAccessRecordConst( allCOUPs[idcoup], ievt0 ); // dependent couplings, vary event-by-event
      for( size_t iicoup = 0; iicoup < nicoup; iicoup++ )
        COUPs[ndcoup + iicoup] = allCOUPs[ndcoup + iicoup]; // independent couplings, fixed for all events
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
#endif
#endif

      // Reset color flows (reset jamp_sv) at the beginning of a new event or event page
      for( int i = 0; i < ncolor; i++ ) { jamp_sv[i] = cxzero_sv(); }

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      // Numerators and denominators for the current event (CUDA) or SIMD event page (C++)
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
#endif

      // *** DIAGRAM 1 OF 123 ***

      // Wavefunction(s) for diagram number 1
      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );

      oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );

      ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][5], +1, w_fp[5], 5 );

      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[0], 0., 0., w_fp[6] );
      FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 0., 0., w_fp[7] );

      // Amplitude(s) for diagram number 1
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[4], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 2 OF 123 ***

      // Wavefunction(s) for diagram number 2
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[4], COUPs[0], 0., 0., w_fp[8] );

      // Amplitude(s) for diagram number 2
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[8], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 3 OF 123 ***

      // Wavefunction(s) for diagram number 3
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], COUPs[0], 0., 0., w_fp[9] );

      // Amplitude(s) for diagram number 3
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[9], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 4 OF 123 ***

      // Wavefunction(s) for diagram number 4
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[5], COUPs[0], 0., 0., w_fp[10] );

      // Amplitude(s) for diagram number 4
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[7], w_fp[10], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 5 OF 123 ***

      // Wavefunction(s) for diagram number 5
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[11] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[6], COUPs[1], cIPD[0], cIPD[1], w_fp[12] );

      // Amplitude(s) for diagram number 5
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 6 OF 123 ***

      // Wavefunction(s) for diagram number 6
      // (none)

      // Amplitude(s) for diagram number 6
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[9], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[12] += amp_sv[0];
      jamp_sv[14] -= amp_sv[0];
      jamp_sv[16] -= amp_sv[0];
      jamp_sv[17] += amp_sv[0];

      // *** DIAGRAM 7 OF 123 ***

      // Wavefunction(s) for diagram number 7
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], COUPs[1], cIPD[0], cIPD[1], w_fp[13] );

      // Amplitude(s) for diagram number 7
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[11], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 7 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 8 OF 123 ***

      // Wavefunction(s) for diagram number 8
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[5], COUPs[1], cIPD[0], cIPD[1], w_fp[14] );

      // Amplitude(s) for diagram number 8
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 8 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 9 OF 123 ***

      // Wavefunction(s) for diagram number 9
      // (none)

      // Amplitude(s) for diagram number 9
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[8], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 9 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[18] += amp_sv[0];
      jamp_sv[20] -= amp_sv[0];
      jamp_sv[22] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];

      // *** DIAGRAM 10 OF 123 ***

      // Wavefunction(s) for diagram number 10
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[15] );

      // Amplitude(s) for diagram number 10
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[14], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 10 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 11 OF 123 ***

      // Wavefunction(s) for diagram number 11
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[6], COUPs[1], cIPD[0], cIPD[1], w_fp[16] );

      // Amplitude(s) for diagram number 11
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[16], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 11 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 12 OF 123 ***

      // Wavefunction(s) for diagram number 12
      // (none)

      // Amplitude(s) for diagram number 12
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[9], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 12 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += amp_sv[0];
      jamp_sv[7] -= amp_sv[0];
      jamp_sv[18] -= amp_sv[0];
      jamp_sv[20] += amp_sv[0];

      // *** DIAGRAM 13 OF 123 ***

      // Wavefunction(s) for diagram number 13
      // (none)

      // Amplitude(s) for diagram number 13
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[16], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 13 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 14 OF 123 ***

      // Wavefunction(s) for diagram number 14
      // (none)

      // Amplitude(s) for diagram number 14
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[8], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 14 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[6] -= amp_sv[0];
      jamp_sv[12] -= amp_sv[0];
      jamp_sv[14] += amp_sv[0];

      // *** DIAGRAM 15 OF 123 ***

      // Wavefunction(s) for diagram number 15
      // (none)

      // Amplitude(s) for diagram number 15
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[16], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 15 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[6] -= amp_sv[0];
      jamp_sv[7] += amp_sv[0];

      // *** DIAGRAM 16 OF 123 ***

      // Wavefunction(s) for diagram number 16
      // (none)

      // Amplitude(s) for diagram number 16
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 16 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[16] += amp_sv[0];
      jamp_sv[17] -= amp_sv[0];
      jamp_sv[22] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];

      // *** DIAGRAM 17 OF 123 ***

      // Wavefunction(s) for diagram number 17
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[12] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[16] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[8] );

      // Amplitude(s) for diagram number 17
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[8], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 17 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] -= amp_sv[0];

      // *** DIAGRAM 18 OF 123 ***

      // Wavefunction(s) for diagram number 18
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], cIPD[0], cIPD[1], w_fp[9] );

      // Amplitude(s) for diagram number 18
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[9], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 18 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[5] -= amp_sv[0];

      // *** DIAGRAM 19 OF 123 ***

      // Wavefunction(s) for diagram number 19
      // (none)

      // Amplitude(s) for diagram number 19
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[12], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 19 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 20 OF 123 ***

      // Wavefunction(s) for diagram number 20
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[0], 0., 0., w_fp[6] );
      FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], 0., 0., w_fp[17] );

      // Amplitude(s) for diagram number 20
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[17], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 20 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[4] -= amp_sv[0];
      jamp_sv[5] += amp_sv[0];

      // *** DIAGRAM 21 OF 123 ***

      // Wavefunction(s) for diagram number 21
      // (none)

      // Amplitude(s) for diagram number 21
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 21 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 22 OF 123 ***

      // Wavefunction(s) for diagram number 22
      // (none)

      // Amplitude(s) for diagram number 22
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[12], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 22 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 23 OF 123 ***

      // Wavefunction(s) for diagram number 23
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], COUPs[0], 0., 0., w_fp[18] );

      // Amplitude(s) for diagram number 23
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[17], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 23 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += amp_sv[0];
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[3] += amp_sv[0];
      jamp_sv[4] -= amp_sv[0];

      // *** DIAGRAM 24 OF 123 ***

      // Wavefunction(s) for diagram number 24
      // (none)

      // Amplitude(s) for diagram number 24
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[8], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 24 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 25 OF 123 ***

      // Wavefunction(s) for diagram number 25
      // (none)

      // Amplitude(s) for diagram number 25
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[12], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 25 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 26 OF 123 ***

      // Wavefunction(s) for diagram number 26
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[19] );

      // Amplitude(s) for diagram number 26
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[19], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 26 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 27 OF 123 ***

      // Wavefunction(s) for diagram number 27
      // (none)

      // Amplitude(s) for diagram number 27
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 27 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] -= amp_sv[0];

      // *** DIAGRAM 28 OF 123 ***

      // Wavefunction(s) for diagram number 28
      // (none)

      // Amplitude(s) for diagram number 28
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[19], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 28 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 29 OF 123 ***

      // Wavefunction(s) for diagram number 29
      // (none)

      // Amplitude(s) for diagram number 29
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[8], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 29 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] -= amp_sv[0];

      // *** DIAGRAM 30 OF 123 ***

      // Wavefunction(s) for diagram number 30
      // (none)

      // Amplitude(s) for diagram number 30
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[19], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 30 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 31 OF 123 ***

      // Wavefunction(s) for diagram number 31
      // (none)

      // Amplitude(s) for diagram number 31
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[17], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 31 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[5] += amp_sv[0];

      // *** DIAGRAM 32 OF 123 ***

      // Wavefunction(s) for diagram number 32
      VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[17] );
      VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[19] );
      VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[8] );

      // Amplitude(s) for diagram number 32
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[17], COUPs[1], &amp_fp[0] );
      jamp_sv[0] += amp_sv[0];
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[5] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[19], COUPs[1], &amp_fp[0] );
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[2] += amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[4] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[8], COUPs[1], &amp_fp[0] );
      jamp_sv[0] -= amp_sv[0];
      jamp_sv[2] += amp_sv[0];
      jamp_sv[4] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];

      // *** DIAGRAM 33 OF 123 ***

      // Wavefunction(s) for diagram number 33
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[12] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[9] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[20] );

      // Amplitude(s) for diagram number 33
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[9], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 33 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[11] -= amp_sv[0];

      // *** DIAGRAM 34 OF 123 ***

      // Wavefunction(s) for diagram number 34
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[1], cIPD[0], cIPD[1], w_fp[21] );

      // Amplitude(s) for diagram number 34
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 34 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[9] -= amp_sv[0];

      // *** DIAGRAM 35 OF 123 ***

      // Wavefunction(s) for diagram number 35
      // (none)

      // Amplitude(s) for diagram number 35
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 35 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 36 OF 123 ***

      // Wavefunction(s) for diagram number 36
      FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], COUPs[1], 0., 0., w_fp[22] );

      // Amplitude(s) for diagram number 36
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[22], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 36 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[9] += amp_sv[0];
      jamp_sv[15] -= amp_sv[0];
      jamp_sv[21] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];

      // *** DIAGRAM 37 OF 123 ***

      // Wavefunction(s) for diagram number 37
      // (none)

      // Amplitude(s) for diagram number 37
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 37 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 38 OF 123 ***

      // Wavefunction(s) for diagram number 38
      // (none)

      // Amplitude(s) for diagram number 38
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 38 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 39 OF 123 ***

      // Wavefunction(s) for diagram number 39
      // (none)

      // Amplitude(s) for diagram number 39
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[22], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 39 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[11] += amp_sv[0];
      jamp_sv[15] -= amp_sv[0];
      jamp_sv[17] += amp_sv[0];
      jamp_sv[21] -= amp_sv[0];

      // *** DIAGRAM 40 OF 123 ***

      // Wavefunction(s) for diagram number 40
      // (none)

      // Amplitude(s) for diagram number 40
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[2], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 40 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 41 OF 123 ***

      // Wavefunction(s) for diagram number 41
      // (none)

      // Amplitude(s) for diagram number 41
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[11], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 41 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 42 OF 123 ***

      // Wavefunction(s) for diagram number 42
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[23] );

      // Amplitude(s) for diagram number 42
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[11], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 42 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[17] -= amp_sv[0];

      // *** DIAGRAM 43 OF 123 ***

      // Wavefunction(s) for diagram number 43
      // (none)

      // Amplitude(s) for diagram number 43
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 43 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[15] -= amp_sv[0];

      // *** DIAGRAM 44 OF 123 ***

      // Wavefunction(s) for diagram number 44
      // (none)

      // Amplitude(s) for diagram number 44
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[14], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 44 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[23] -= amp_sv[0];

      // *** DIAGRAM 45 OF 123 ***

      // Wavefunction(s) for diagram number 45
      // (none)

      // Amplitude(s) for diagram number 45
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[14], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 45 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[21] -= amp_sv[0];

      // *** DIAGRAM 46 OF 123 ***

      // Wavefunction(s) for diagram number 46
      // (none)

      // Amplitude(s) for diagram number 46
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[23], w_fp[2], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 46 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 47 OF 123 ***

      // Wavefunction(s) for diagram number 47
      // (none)

      // Amplitude(s) for diagram number 47
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[22], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 47 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[9] += amp_sv[0];
      jamp_sv[11] -= amp_sv[0];
      jamp_sv[17] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];

      // *** DIAGRAM 48 OF 123 ***

      // Wavefunction(s) for diagram number 48
      // (none)

      // Amplitude(s) for diagram number 48
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[17], COUPs[1], &amp_fp[0] );
      jamp_sv[9] += amp_sv[0];
      jamp_sv[11] -= amp_sv[0];
      jamp_sv[17] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[19], COUPs[1], &amp_fp[0] );
      jamp_sv[11] -= amp_sv[0];
      jamp_sv[15] += amp_sv[0];
      jamp_sv[17] -= amp_sv[0];
      jamp_sv[21] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[8], COUPs[1], &amp_fp[0] );
      jamp_sv[9] -= amp_sv[0];
      jamp_sv[15] += amp_sv[0];
      jamp_sv[21] += amp_sv[0];
      jamp_sv[23] -= amp_sv[0];

      // *** DIAGRAM 49 OF 123 ***

      // Wavefunction(s) for diagram number 49
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[0], 0., 0., w_fp[12] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], COUPs[1], cIPD[0], cIPD[1], w_fp[22] );

      // Amplitude(s) for diagram number 49
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[9], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 49 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 50 OF 123 ***

      // Wavefunction(s) for diagram number 50
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[5], COUPs[0], 0., 0., w_fp[23] );

      // Amplitude(s) for diagram number 50
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 50 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[6] += amp_sv[0];
      jamp_sv[8] -= amp_sv[0];
      jamp_sv[10] -= amp_sv[0];
      jamp_sv[11] += amp_sv[0];

      // *** DIAGRAM 51 OF 123 ***

      // Wavefunction(s) for diagram number 51
      // (none)

      // Amplitude(s) for diagram number 51
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[9], w_fp[12], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 51 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 52 OF 123 ***

      // Wavefunction(s) for diagram number 52
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[12], COUPs[1], cIPD[0], cIPD[1], w_fp[20] );

      // Amplitude(s) for diagram number 52
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 52 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 53 OF 123 ***

      // Wavefunction(s) for diagram number 53
      // (none)

      // Amplitude(s) for diagram number 53
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 53 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] += amp_sv[0];
      jamp_sv[13] -= amp_sv[0];
      jamp_sv[19] -= amp_sv[0];
      jamp_sv[22] += amp_sv[0];

      // *** DIAGRAM 54 OF 123 ***

      // Wavefunction(s) for diagram number 54
      // (none)

      // Amplitude(s) for diagram number 54
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[14], w_fp[12], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 54 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 55 OF 123 ***

      // Wavefunction(s) for diagram number 55
      // (none)

      // Amplitude(s) for diagram number 55
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 55 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[12] -= amp_sv[0];
      jamp_sv[13] += amp_sv[0];

      // *** DIAGRAM 56 OF 123 ***

      // Wavefunction(s) for diagram number 56
      // (none)

      // Amplitude(s) for diagram number 56
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 56 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[10] += amp_sv[0];
      jamp_sv[11] -= amp_sv[0];
      jamp_sv[20] -= amp_sv[0];
      jamp_sv[21] += amp_sv[0];

      // *** DIAGRAM 57 OF 123 ***

      // Wavefunction(s) for diagram number 57
      // (none)

      // Amplitude(s) for diagram number 57
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[18], w_fp[7], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 57 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 58 OF 123 ***

      // Wavefunction(s) for diagram number 58
      // (none)

      // Amplitude(s) for diagram number 58
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 59 OF 123 ***

      // Wavefunction(s) for diagram number 59
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[12], w_fp[1], COUPs[0], 0., 0., w_fp[21] );

      // Amplitude(s) for diagram number 59
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[21], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 59 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 60 OF 123 ***

      // Wavefunction(s) for diagram number 60
      // (none)

      // Amplitude(s) for diagram number 60
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[23], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 60 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 61 OF 123 ***

      // Wavefunction(s) for diagram number 61
      // (none)

      // Amplitude(s) for diagram number 61
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[21], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 61 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[19] += amp_sv[0];
      jamp_sv[20] -= amp_sv[0];
      jamp_sv[21] += amp_sv[0];
      jamp_sv[22] -= amp_sv[0];

      // *** DIAGRAM 62 OF 123 ***

      // Wavefunction(s) for diagram number 62
      // (none)

      // Amplitude(s) for diagram number 62
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 62 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 63 OF 123 ***

      // Wavefunction(s) for diagram number 63
      // (none)

      // Amplitude(s) for diagram number 63
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[21], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 63 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += amp_sv[0];
      jamp_sv[6] -= amp_sv[0];
      jamp_sv[8] += amp_sv[0];
      jamp_sv[12] -= amp_sv[0];

      // *** DIAGRAM 64 OF 123 ***

      // Wavefunction(s) for diagram number 64
      // (none)

      // Amplitude(s) for diagram number 64
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[20], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 64 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 65 OF 123 ***

      // Wavefunction(s) for diagram number 65
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[5], COUPs[0], 0., 0., w_fp[20] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], COUPs[1], cIPD[0], cIPD[1], w_fp[21] );

      // Amplitude(s) for diagram number 65
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 65 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 66 OF 123 ***

      // Wavefunction(s) for diagram number 66
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[4], COUPs[0], 0., 0., w_fp[22] );

      // Amplitude(s) for diagram number 66
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[22], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 66 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[7] += amp_sv[0];
      jamp_sv[8] -= amp_sv[0];
      jamp_sv[9] += amp_sv[0];
      jamp_sv[10] -= amp_sv[0];

      // *** DIAGRAM 67 OF 123 ***

      // Wavefunction(s) for diagram number 67
      // (none)

      // Amplitude(s) for diagram number 67
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[9], w_fp[20], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 67 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 68 OF 123 ***

      // Wavefunction(s) for diagram number 68
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[20], COUPs[1], cIPD[0], cIPD[1], w_fp[23] );

      // Amplitude(s) for diagram number 68
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[23], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 68 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 69 OF 123 ***

      // Wavefunction(s) for diagram number 69
      // (none)

      // Amplitude(s) for diagram number 69
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[22], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 69 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[5] += amp_sv[0];
      jamp_sv[13] -= amp_sv[0];
      jamp_sv[16] += amp_sv[0];
      jamp_sv[19] -= amp_sv[0];

      // *** DIAGRAM 70 OF 123 ***

      // Wavefunction(s) for diagram number 70
      // (none)

      // Amplitude(s) for diagram number 70
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[11], w_fp[20], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 70 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 71 OF 123 ***

      // Wavefunction(s) for diagram number 71
      // (none)

      // Amplitude(s) for diagram number 71
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 71 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];
      jamp_sv[18] -= amp_sv[0];
      jamp_sv[19] += amp_sv[0];

      // *** DIAGRAM 72 OF 123 ***

      // Wavefunction(s) for diagram number 72
      // (none)

      // Amplitude(s) for diagram number 72
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 72 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[8] += amp_sv[0];
      jamp_sv[9] -= amp_sv[0];
      jamp_sv[14] -= amp_sv[0];
      jamp_sv[15] += amp_sv[0];

      // *** DIAGRAM 73 OF 123 ***

      // Wavefunction(s) for diagram number 73
      // (none)

      // Amplitude(s) for diagram number 73
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[6], w_fp[7], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 73 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 74 OF 123 ***

      // Wavefunction(s) for diagram number 74
      // (none)

      // Amplitude(s) for diagram number 74
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 75 OF 123 ***

      // Wavefunction(s) for diagram number 75
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], COUPs[0], 0., 0., w_fp[12] );

      // Amplitude(s) for diagram number 75
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[12], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 75 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 76 OF 123 ***

      // Wavefunction(s) for diagram number 76
      // (none)

      // Amplitude(s) for diagram number 76
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[22], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 76 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 77 OF 123 ***

      // Wavefunction(s) for diagram number 77
      // (none)

      // Amplitude(s) for diagram number 77
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[12], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 77 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[13] += amp_sv[0];
      jamp_sv[14] -= amp_sv[0];
      jamp_sv[15] += amp_sv[0];
      jamp_sv[16] -= amp_sv[0];

      // *** DIAGRAM 78 OF 123 ***

      // Wavefunction(s) for diagram number 78
      // (none)

      // Amplitude(s) for diagram number 78
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 78 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 79 OF 123 ***

      // Wavefunction(s) for diagram number 79
      // (none)

      // Amplitude(s) for diagram number 79
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[12], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 79 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += amp_sv[0];
      jamp_sv[7] -= amp_sv[0];
      jamp_sv[10] += amp_sv[0];
      jamp_sv[18] -= amp_sv[0];

      // *** DIAGRAM 80 OF 123 ***

      // Wavefunction(s) for diagram number 80
      // (none)

      // Amplitude(s) for diagram number 80
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 80 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 81 OF 123 ***

      // Wavefunction(s) for diagram number 81
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[9], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[23] );

      // Amplitude(s) for diagram number 81
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[23], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 81 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[7] -= amp_sv[0];

      // *** DIAGRAM 82 OF 123 ***

      // Wavefunction(s) for diagram number 82
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[15], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[12] );

      // Amplitude(s) for diagram number 82
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 82 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[10] -= amp_sv[0];

      // *** DIAGRAM 83 OF 123 ***

      // Wavefunction(s) for diagram number 83
      // (none)

      // Amplitude(s) for diagram number 83
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[23], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 83 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[6] -= amp_sv[0];

      // *** DIAGRAM 84 OF 123 ***

      // Wavefunction(s) for diagram number 84
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[13], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[21] );

      // Amplitude(s) for diagram number 84
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 84 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[8] -= amp_sv[0];

      // *** DIAGRAM 85 OF 123 ***

      // Wavefunction(s) for diagram number 85
      // (none)

      // Amplitude(s) for diagram number 85
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 85 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 86 OF 123 ***

      // Wavefunction(s) for diagram number 86
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], COUPs[0], 0., 0., w_fp[23] );

      // Amplitude(s) for diagram number 86
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[23], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 86 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[6] += amp_sv[0];
      jamp_sv[7] -= amp_sv[0];
      jamp_sv[9] -= amp_sv[0];
      jamp_sv[11] += amp_sv[0];

      // *** DIAGRAM 87 OF 123 ***

      // Wavefunction(s) for diagram number 87
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[16], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[22] );

      // Amplitude(s) for diagram number 87
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[11], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 87 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[16] -= amp_sv[0];

      // *** DIAGRAM 88 OF 123 ***

      // Wavefunction(s) for diagram number 88
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[11], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[20] );

      // Amplitude(s) for diagram number 88
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[20], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 88 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[13] -= amp_sv[0];

      // *** DIAGRAM 89 OF 123 ***

      // Wavefunction(s) for diagram number 89
      // (none)

      // Amplitude(s) for diagram number 89
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[14], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 89 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[22] -= amp_sv[0];

      // *** DIAGRAM 90 OF 123 ***

      // Wavefunction(s) for diagram number 90
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[14], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[24] );

      // Amplitude(s) for diagram number 90
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[24], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 90 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[19] -= amp_sv[0];

      // *** DIAGRAM 91 OF 123 ***

      // Wavefunction(s) for diagram number 91
      // (none)

      // Amplitude(s) for diagram number 91
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 91 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 92 OF 123 ***

      // Wavefunction(s) for diagram number 92
      // (none)

      // Amplitude(s) for diagram number 92
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[23], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 92 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];
      jamp_sv[16] -= amp_sv[0];
      jamp_sv[22] += amp_sv[0];

      // *** DIAGRAM 93 OF 123 ***

      // Wavefunction(s) for diagram number 93
      // (none)

      // Amplitude(s) for diagram number 93
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], w_fp[7], w_fp[5], COUPs[2], &amp_fp[0] );
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 94 OF 123 ***

      // Wavefunction(s) for diagram number 94
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[6], COUPs[0], 0., 0., w_fp[22] );

      // Amplitude(s) for diagram number 94
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[5], w_fp[22], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 94 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 95 OF 123 ***

      // Wavefunction(s) for diagram number 95
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[7], COUPs[0], 0., 0., w_fp[25] );

      // Amplitude(s) for diagram number 95
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[6], w_fp[5], w_fp[25], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 95 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 96 OF 123 ***

      // Wavefunction(s) for diagram number 96
      // (none)

      // Amplitude(s) for diagram number 96
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[22], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 96 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[18] += amp_sv[0];
      jamp_sv[19] -= amp_sv[0];
      jamp_sv[21] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];

      // *** DIAGRAM 97 OF 123 ***

      // Wavefunction(s) for diagram number 97
      // (none)

      // Amplitude(s) for diagram number 97
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[24], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 97 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 98 OF 123 ***

      // Wavefunction(s) for diagram number 98
      // (none)

      // Amplitude(s) for diagram number 98
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[22], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 98 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[8] -= amp_sv[0];
      jamp_sv[14] += amp_sv[0];

      // *** DIAGRAM 99 OF 123 ***

      // Wavefunction(s) for diagram number 99
      // (none)

      // Amplitude(s) for diagram number 99
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 99 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 100 OF 123 ***

      // Wavefunction(s) for diagram number 100
      // (none)

      // Amplitude(s) for diagram number 100
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], w_fp[7], w_fp[4], COUPs[2], &amp_fp[0] );
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 101 OF 123 ***

      // Wavefunction(s) for diagram number 101
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[18], COUPs[0], 0., 0., w_fp[6] );

      // Amplitude(s) for diagram number 101
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[4], w_fp[6], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 101 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 102 OF 123 ***

      // Wavefunction(s) for diagram number 102
      // (none)

      // Amplitude(s) for diagram number 102
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[4], w_fp[25], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 102 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 103 OF 123 ***

      // Wavefunction(s) for diagram number 103
      // (none)

      // Amplitude(s) for diagram number 103
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 103 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[12] += amp_sv[0];
      jamp_sv[13] -= amp_sv[0];
      jamp_sv[15] -= amp_sv[0];
      jamp_sv[17] += amp_sv[0];

      // *** DIAGRAM 104 OF 123 ***

      // Wavefunction(s) for diagram number 104
      // (none)

      // Amplitude(s) for diagram number 104
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 104 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 105 OF 123 ***

      // Wavefunction(s) for diagram number 105
      // (none)

      // Amplitude(s) for diagram number 105
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[6], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 105 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += amp_sv[0];
      jamp_sv[4] -= amp_sv[0];
      jamp_sv[10] -= amp_sv[0];
      jamp_sv[20] += amp_sv[0];

      // *** DIAGRAM 106 OF 123 ***

      // Wavefunction(s) for diagram number 106
      // (none)

      // Amplitude(s) for diagram number 106
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[18], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 106 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 107 OF 123 ***

      // Wavefunction(s) for diagram number 107
      // (none)

      // Amplitude(s) for diagram number 107
      VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[7], w_fp[10], COUPs[2], &amp_fp[0] );
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 108 OF 123 ***

      // Wavefunction(s) for diagram number 108
      // (none)

      // Amplitude(s) for diagram number 108
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[10], w_fp[25], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 108 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 109 OF 123 ***

      // Wavefunction(s) for diagram number 109
      // (none)

      // Amplitude(s) for diagram number 109
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[7], w_fp[23], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 109 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 110 OF 123 ***

      // Wavefunction(s) for diagram number 110
      // (none)

      // Amplitude(s) for diagram number 110
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[20], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 110 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[12] -= amp_sv[0];

      // *** DIAGRAM 111 OF 123 ***

      // Wavefunction(s) for diagram number 111
      // (none)

      // Amplitude(s) for diagram number 111
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[11], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 111 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[14] -= amp_sv[0];

      // *** DIAGRAM 112 OF 123 ***

      // Wavefunction(s) for diagram number 112
      // (none)

      // Amplitude(s) for diagram number 112
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[24], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 112 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[18] -= amp_sv[0];

      // *** DIAGRAM 113 OF 123 ***

      // Wavefunction(s) for diagram number 113
      // (none)

      // Amplitude(s) for diagram number 113
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[14], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 113 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[20] -= amp_sv[0];

      // *** DIAGRAM 114 OF 123 ***

      // Wavefunction(s) for diagram number 114
      VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[12] );
      VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[24] );
      VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[21] );

      // Amplitude(s) for diagram number 114
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[7], w_fp[5], COUPs[0], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[7], w_fp[5], COUPs[0], &amp_fp[0] );
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[7], w_fp[5], COUPs[0], &amp_fp[0] );
      jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 115 OF 123 ***

      // Wavefunction(s) for diagram number 115
      // (none)

      // Amplitude(s) for diagram number 115
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[12], COUPs[1], &amp_fp[0] );
      jamp_sv[18] += amp_sv[0];
      jamp_sv[19] -= amp_sv[0];
      jamp_sv[21] -= amp_sv[0];
      jamp_sv[23] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[19] -= amp_sv[0];
      jamp_sv[20] += amp_sv[0];
      jamp_sv[21] -= amp_sv[0];
      jamp_sv[22] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[14], w_fp[21], COUPs[1], &amp_fp[0] );
      jamp_sv[18] -= amp_sv[0];
      jamp_sv[20] += amp_sv[0];
      jamp_sv[22] += amp_sv[0];
      jamp_sv[23] -= amp_sv[0];

      // *** DIAGRAM 116 OF 123 ***

      // Wavefunction(s) for diagram number 116
      // (none)

      // Amplitude(s) for diagram number 116
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[12], COUPs[1], &amp_fp[0] );
      jamp_sv[0] += amp_sv[0];
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[8] -= amp_sv[0];
      jamp_sv[14] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[6] += amp_sv[0];
      jamp_sv[8] -= amp_sv[0];
      jamp_sv[12] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[2], w_fp[21], COUPs[1], &amp_fp[0] );
      jamp_sv[0] -= amp_sv[0];
      jamp_sv[6] += amp_sv[0];
      jamp_sv[12] += amp_sv[0];
      jamp_sv[14] -= amp_sv[0];

      // *** DIAGRAM 117 OF 123 ***

      // Wavefunction(s) for diagram number 117
      VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 0., 0., w_fp[21] );
      VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 0., 0., w_fp[13] );
      VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 0., 0., w_fp[24] );

      // Amplitude(s) for diagram number 117
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[7], w_fp[4], COUPs[0], &amp_fp[0] );
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[7], w_fp[4], COUPs[0], &amp_fp[0] );
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[7], w_fp[4], COUPs[0], &amp_fp[0] );
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 118 OF 123 ***

      // Wavefunction(s) for diagram number 118
      // (none)

      // Amplitude(s) for diagram number 118
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[21], COUPs[1], &amp_fp[0] );
      jamp_sv[12] += amp_sv[0];
      jamp_sv[13] -= amp_sv[0];
      jamp_sv[15] -= amp_sv[0];
      jamp_sv[17] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[13], COUPs[1], &amp_fp[0] );
      jamp_sv[13] -= amp_sv[0];
      jamp_sv[14] += amp_sv[0];
      jamp_sv[15] -= amp_sv[0];
      jamp_sv[16] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[12] -= amp_sv[0];
      jamp_sv[14] += amp_sv[0];
      jamp_sv[16] += amp_sv[0];
      jamp_sv[17] -= amp_sv[0];

      // *** DIAGRAM 119 OF 123 ***

      // Wavefunction(s) for diagram number 119
      // (none)

      // Amplitude(s) for diagram number 119
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[21], COUPs[1], &amp_fp[0] );
      jamp_sv[1] += amp_sv[0];
      jamp_sv[4] -= amp_sv[0];
      jamp_sv[10] -= amp_sv[0];
      jamp_sv[20] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[13], COUPs[1], &amp_fp[0] );
      jamp_sv[4] -= amp_sv[0];
      jamp_sv[7] += amp_sv[0];
      jamp_sv[10] -= amp_sv[0];
      jamp_sv[18] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[2], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[7] += amp_sv[0];
      jamp_sv[18] += amp_sv[0];
      jamp_sv[20] -= amp_sv[0];

      // *** DIAGRAM 120 OF 123 ***

      // Wavefunction(s) for diagram number 120
      VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[24] );
      VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[15] );
      VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 0., 0., w_fp[13] );

      // Amplitude(s) for diagram number 120
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[6] += amp_sv[0];
      jamp_sv[7] -= amp_sv[0];
      jamp_sv[9] -= amp_sv[0];
      jamp_sv[11] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[15], COUPs[1], &amp_fp[0] );
      jamp_sv[7] -= amp_sv[0];
      jamp_sv[8] += amp_sv[0];
      jamp_sv[9] -= amp_sv[0];
      jamp_sv[10] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[13], COUPs[1], &amp_fp[0] );
      jamp_sv[6] -= amp_sv[0];
      jamp_sv[8] += amp_sv[0];
      jamp_sv[10] += amp_sv[0];
      jamp_sv[11] -= amp_sv[0];

      // *** DIAGRAM 121 OF 123 ***

      // Wavefunction(s) for diagram number 121
      // (none)

      // Amplitude(s) for diagram number 121
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[24], COUPs[1], &amp_fp[0] );
      jamp_sv[3] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];
      jamp_sv[16] -= amp_sv[0];
      jamp_sv[22] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[15], COUPs[1], &amp_fp[0] );
      jamp_sv[5] -= amp_sv[0];
      jamp_sv[13] += amp_sv[0];
      jamp_sv[16] -= amp_sv[0];
      jamp_sv[19] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[13], COUPs[1], &amp_fp[0] );
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[13] += amp_sv[0];
      jamp_sv[19] += amp_sv[0];
      jamp_sv[22] -= amp_sv[0];

      // *** DIAGRAM 122 OF 123 ***

      // Wavefunction(s) for diagram number 122
      // (none)

      // Amplitude(s) for diagram number 122
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[15], w_fp[1], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[13], w_fp[1], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 123 OF 123 ***

      // Wavefunction(s) for diagram number 123
      // (none)

      // Amplitude(s) for diagram number 123
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[17], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[19], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[7], COUPs[0], &amp_fp[0] );
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called CPPProcess::matrix_1_gg_ttxgg()?)

      // The color denominators (initialize all array elements, with ncolor=24)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype denom[ncolor] = { 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54 }; // 1-D array[24]

      // The color matrix (initialize all array elements, with ncolor=24)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype cf[ncolor][ncolor] = {
        { 512, -64, -64, 8, 8, 80, -64, 8, 8, -1, -1, -10, 8, -1, 80, -10, 71, 62, -1, -10, -10, 62, 62, -28 },
        { -64, 512, 8, 80, -64, 8, 8, -64, -1, -10, 8, -1, -1, -10, -10, 62, 62, -28, 8, -1, 80, -10, 71, 62 },
        { -64, 8, 512, -64, 80, 8, 8, -1, 80, -10, 71, 62, -64, 8, 8, -1, -1, -10, -10, -1, 62, -28, -10, 62 },
        { 8, 80, -64, 512, 8, -64, -1, -10, -10, 62, 62, -28, 8, -64, -1, -10, 8, -1, -1, 8, 71, 62, 80, -10 },
        { 8, -64, 80, 8, 512, -64, -1, 8, 71, 62, 80, -10, -10, -1, 62, -28, -10, 62, -64, 8, 8, -1, -1, -10 },
        { 80, 8, 8, -64, -64, 512, -10, -1, 62, -28, -10, 62, -1, 8, 71, 62, 80, -10, 8, -64, -1, -10, 8, -1 },
        { -64, 8, 8, -1, -1, -10, 512, -64, -64, 8, 8, 80, 80, -10, 8, -1, 62, 71, -10, 62, -1, -10, -28, 62 },
        { 8, -64, -1, -10, 8, -1, -64, 512, 8, 80, -64, 8, -10, 62, -1, -10, -28, 62, 80, -10, 8, -1, 62, 71 },
        { 8, -1, 80, -10, 71, 62, -64, 8, 512, -64, 80, 8, 8, -1, -64, 8, -10, -1, 62, -28, -10, -1, 62, -10 },
        { -1, -10, -10, 62, 62, -28, 8, 80, -64, 512, 8, -64, -1, -10, 8, -64, -1, 8, 71, 62, -1, 8, -10, 80 },
        { -1, 8, 71, 62, 80, -10, 8, -64, 80, 8, 512, -64, 62, -28, -10, -1, 62, -10, 8, -1, -64, 8, -10, -1 },
        { -10, -1, 62, -28, -10, 62, 80, 8, 8, -64, -64, 512, 71, 62, -1, 8, -10, 80, -1, -10, 8, -64, -1, 8 },
        { 8, -1, -64, 8, -10, -1, 80, -10, 8, -1, 62, 71, 512, -64, -64, 8, 8, 80, 62, -10, -28, 62, -1, -10 },
        { -1, -10, 8, -64, -1, 8, -10, 62, -1, -10, -28, 62, -64, 512, 8, 80, -64, 8, -10, 80, 62, 71, 8, -1 },
        { 80, -10, 8, -1, 62, 71, 8, -1, -64, 8, -10, -1, -64, 8, 512, -64, 80, 8, -28, 62, 62, -10, -10, -1 },
        { -10, 62, -1, -10, -28, 62, -1, -10, 8, -64, -1, 8, 8, 80, -64, 512, 8, -64, 62, 71, -10, 80, -1, 8 },
        { 71, 62, -1, 8, -10, 80, 62, -28, -10, -1, 62, -10, 8, -64, 80, 8, 512, -64, -1, 8, -10, -1, -64, 8 },
        { 62, -28, -10, -1, 62, -10, 71, 62, -1, 8, -10, 80, 80, 8, 8, -64, -64, 512, -10, -1, -1, 8, 8, -64 },
        { -1, 8, -10, -1, -64, 8, -10, 80, 62, 71, 8, -1, 62, -10, -28, 62, -1, -10, 512, -64, -64, 8, 8, 80 },
        { -10, -1, -1, 8, 8, -64, 62, -10, -28, 62, -1, -10, -10, 80, 62, 71, 8, -1, -64, 512, 8, 80, -64, 8 },
        { -10, 80, 62, 71, 8, -1, -1, 8, -10, -1, -64, 8, -28, 62, 62, -10, -10, -1, -64, 8, 512, -64, 80, 8 },
        { 62, -10, -28, 62, -1, -10, -10, -1, -1, 8, 8, -64, 62, 71, -10, 80, -1, 8, 8, 80, -64, 512, 8, -64 },
        { 62, 71, -10, 80, -1, 8, -28, 62, 62, -10, -10, -1, -1, 8, -10, -1, -64, 8, 8, -64, 80, 8, 512, -64 },
        { -28, 62, 62, -10, -10, -1, 62, 71, -10, 80, -1, 8, -10, -1, -1, 8, 8, -64, 80, 8, 8, -64, -64, 512 } }; // 2-D array[24][24]

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
      // Use the property that M is a real matrix (see #475):
      // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
      for( int icol = 0; icol < ncolor; icol++ )
      {
        // In addition, for C++ use the property that M is symmetric (see #475):
        // we gain (a factor 2?) in speed here as we only loop over the up diagonal part of the matrix!
        // For CUDA, the old implementation is (surprisingly?!) faster expecially for complex physics and large grids
#ifndef __CUDACC__
        // Diagonal terms
        fptype_sv jampRi_sv = cxreal( jamp_sv[icol] );
        fptype_sv jampIi_sv = cximag( jamp_sv[icol] );
        fptype_sv ztempR_sv = cf[icol][icol] * jampRi_sv;
        fptype_sv ztempI_sv = cf[icol][icol] * jampIi_sv;
        // Off-diagonal terms
        for( int jcol = icol + 1; jcol < ncolor; jcol++ )
        {
          fptype_sv jampRj_sv = cxreal( jamp_sv[jcol] );
          fptype_sv jampIj_sv = cximag( jamp_sv[jcol] );
          ztempR_sv += 2 * cf[icol][jcol] * jampRj_sv;
          ztempI_sv += 2 * cf[icol][jcol] * jampIj_sv;
        }
        deltaMEs += ( jampRi_sv * ztempR_sv + jampIi_sv * ztempI_sv ) / denom[icol];
#else
        cxtype_sv ztemp_sv = cxzero_sv();
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp_sv += cf[icol][jcol] * jamp_sv[jcol];
        deltaMEs += ( cxreal( ztemp_sv ) * cxreal( jamp_sv[icol] ) + cximag( ztemp_sv ) * cximag( jamp_sv[icol] ) ) / denom[icol];
#endif
      }

      // *** STORE THE RESULTS ***

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv += deltaMEs; // fix #435
      /*
#ifdef __CUDACC__
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", blockDim.x * blockIdx.x + threadIdx.x, ihel, MEs_sv );
#else
#ifdef MGONGPU_CPPSIMD
      if( cNGoodHel > 0 )
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
          printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", ipagV * neppV + ieppV, ihel, MEs_sv[ieppV] );
#else
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%6d ihel=%2d me_running=%f\n", ipagV, ihel, MEs_sv );
#endif
#endif
      */
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( bool verbose,
                          bool debug )
    : m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_PARAM
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
    static constexpr short tHel[ncomb][mgOnGpu::npar] = {
      { -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, 1 },
      { -1, -1, -1, -1, 1, -1 },
      { -1, -1, -1, -1, 1, 1 },
      { -1, -1, -1, 1, -1, -1 },
      { -1, -1, -1, 1, -1, 1 },
      { -1, -1, -1, 1, 1, -1 },
      { -1, -1, -1, 1, 1, 1 },
      { -1, -1, 1, -1, -1, -1 },
      { -1, -1, 1, -1, -1, 1 },
      { -1, -1, 1, -1, 1, -1 },
      { -1, -1, 1, -1, 1, 1 },
      { -1, -1, 1, 1, -1, -1 },
      { -1, -1, 1, 1, -1, 1 },
      { -1, -1, 1, 1, 1, -1 },
      { -1, -1, 1, 1, 1, 1 },
      { -1, 1, -1, -1, -1, -1 },
      { -1, 1, -1, -1, -1, 1 },
      { -1, 1, -1, -1, 1, -1 },
      { -1, 1, -1, -1, 1, 1 },
      { -1, 1, -1, 1, -1, -1 },
      { -1, 1, -1, 1, -1, 1 },
      { -1, 1, -1, 1, 1, -1 },
      { -1, 1, -1, 1, 1, 1 },
      { -1, 1, 1, -1, -1, -1 },
      { -1, 1, 1, -1, -1, 1 },
      { -1, 1, 1, -1, 1, -1 },
      { -1, 1, 1, -1, 1, 1 },
      { -1, 1, 1, 1, -1, -1 },
      { -1, 1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1, 1 },
      { 1, -1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, -1, 1 },
      { 1, -1, -1, -1, 1, -1 },
      { 1, -1, -1, -1, 1, 1 },
      { 1, -1, -1, 1, -1, -1 },
      { 1, -1, -1, 1, -1, 1 },
      { 1, -1, -1, 1, 1, -1 },
      { 1, -1, -1, 1, 1, 1 },
      { 1, -1, 1, -1, -1, -1 },
      { 1, -1, 1, -1, -1, 1 },
      { 1, -1, 1, -1, 1, -1 },
      { 1, -1, 1, -1, 1, 1 },
      { 1, -1, 1, 1, -1, -1 },
      { 1, -1, 1, 1, -1, 1 },
      { 1, -1, 1, 1, 1, -1 },
      { 1, -1, 1, 1, 1, 1 },
      { 1, 1, -1, -1, -1, -1 },
      { 1, 1, -1, -1, -1, 1 },
      { 1, 1, -1, -1, 1, -1 },
      { 1, 1, -1, -1, 1, 1 },
      { 1, 1, -1, 1, -1, -1 },
      { 1, 1, -1, 1, -1, 1 },
      { 1, 1, -1, 1, 1, -1 },
      { 1, 1, -1, 1, 1, 1 },
      { 1, 1, 1, -1, -1, -1 },
      { 1, 1, 1, -1, -1, 1 },
      { 1, 1, 1, -1, 1, -1 },
      { 1, 1, 1, -1, 1, 1 },
      { 1, 1, 1, 1, -1, -1 },
      { 1, 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1, 1 } };
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * mgOnGpu::npar * sizeof( short ) ) );
#else
    memcpy( cHel, tHel, ncomb * mgOnGpu::npar * sizeof( short ) );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_PARAM
  // Initialize process (with parameters read from user cards)
  void
  CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    //m_pars->setDependentParameters(); // now computed event-by-event (running alphas #373)
    //m_pars->setDependentCouplings(); // now computed event-by-event (running alphas #373)
    if( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      //m_pars->printDependentParameters(); // now computed event-by-event (running alphas #373)
      //m_pars->printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->mdl_MT );
    m_masses.push_back( m_pars->mdl_MT );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    const fptype tIPD[2] = { (fptype)m_pars->mdl_MT, (fptype)m_pars->mdl_WT };
    //const cxtype tIPC[0] = { ... }; // nicoup=0
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof( fptype ) ) );
    //checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 0 * sizeof( cxtype ) ) ); // nicoup=0
#else
    memcpy( cIPD, tIPD, 2 * sizeof( fptype ) );
    //memcpy( cIPC, tIPC, 0 * sizeof( cxtype ) ); // nicoup=0
#endif
    //for ( i=0; i<2; i++ ) std::cout << std::setprecision(17) << "tIPD[i] = " << tIPD[i] << std::endl;
  }
#else
  // Initialize process (with hardcoded parameters)
  void
  CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if( m_verbose )
    {
      Parameters_sm::printIndependentParameters();
      Parameters_sm::printIndependentCouplings();
      //Parameters_sm::printDependentParameters(); // now computed event-by-event (running alphas #373)
      //Parameters_sm::printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::mdl_MT );
    m_masses.push_back( Parameters_sm::mdl_MT );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string
  CPPProcess::getCompiler()
  {
    std::stringstream out;
    // CUDA version (NVCC)
    // [Use __NVCC__ instead of __CUDACC__ here!]
    // [This tests if 'nvcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'nvcc --compiler-options -dM -E dummy.c | grep CUDA': see https://stackoverflow.com/a/53713712]
#ifdef __NVCC__
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
    out << "icx " << __INTEL_LLVM_COMPILER;
#ifdef __NVCC__
    out << ", ";
#else
    out << " (";
#endif
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
#ifdef __APPLE__
    out << "Apple clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype( &pclose )> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
#endif
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
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

  __global__ void /* clang-format off */
  computeDependentCouplings( const fptype* allgs, // input: Gs[nevt]
                             fptype* allcouplings // output: couplings[nevt*ndcoup*2]
#ifndef __CUDACC__
                             , const int nevt     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
  ) /* clang-format on */
  {
#ifdef __CUDACC__
    using namespace mg5amcGpu;
    using G_ACCESS = DeviceAccessGs;
    using C_ACCESS = DeviceAccessCouplings;
    G2COUP<G_ACCESS, C_ACCESS>( allgs, allcouplings );
#else
    using namespace mg5amcCpu;
    using G_ACCESS = HostAccessGs;
    using C_ACCESS = HostAccessCouplings;
    for( int ipagV = 0; ipagV < nevt / neppV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      const fptype* gs = MemoryAccessGs::ieventAccessRecordConst( allgs, ievt0 );
      fptype* couplings = MemoryAccessCouplings::ieventAccessRecord( allcouplings, ievt0 );
      G2COUP<G_ACCESS, C_ACCESS>( gs, couplings );
    }
#endif
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__ /* clang-format off */
  __global__ void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel )           // output: isGoodHel[ncomb] - device array
  { /* clang-format on */
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr unsigned int channelId = 0; // disable single-diagram channel enhancement
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs );
#endif
      if( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                       fptype* allNumerators,      // output: multichannel numerators[nevt], running_sum_over_helicities
                       fptype* allDenominators,    // output: multichannel denominators[nevt], running_sum_over_helicities
#endif
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - device array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    const int maxtry0 = ( neppV > 16 ? neppV : 16 ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    fptype allMEsLast[maxtry0] = { 0 };              // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    const int maxtry = std::min( maxtry0, nevt );    // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for( int ievt = 0; ievt < maxtry; ++ievt )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr unsigned int channelId = 0; // disable single-diagram channel enhancement
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId, maxtry );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, maxtry );
#endif
      for( int ievt = 0; ievt < maxtry; ++ievt )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        const bool differs = ( allMEs[ievt] != allMEsLast[ievt] );
        if( differs )
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

  void
  sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel = 0;           // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
    int goodHel[ncomb] = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if( isGoodHel[ihel] )
      {
        //goodHel[nGoodHel[0]] = ihel; // FIXME: assume process.nprocesses == 1 for the moment
        //nGoodHel[0]++; // FIXME: assume process.nprocesses == 1 for the moment
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof( int ) ) ); // FIXME: assume process.nprocesses == 1 for the moment
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb * sizeof( int ) ) );
#else
    cNGoodHel = nGoodHel;
    for( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  __global__ void /* clang-format off */
  sigmaKin( const fptype* allmomenta,      // input: momenta[nevt*npar*4]
            const fptype* allcouplings,    // input: couplings[nevt*ndcoup*2]
            fptype* allMEs                 // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
            , fptype* allNumerators        // output: multichannel numerators[nevt], running_sum_over_helicities
            , fptype* allDenominators      // output: multichannel denominators[nevt], running_sum_over_helicities
            , const unsigned int channelId // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
#ifndef __CUDACC__
            , const int nevt               // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
            ) /* clang-format on */
  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    constexpr int nprocesses = 1;
    static_assert( nprocesses == 1, "Assume nprocesses == 1" ); // FIXME (#343): assume nprocesses == 1
    constexpr int denominators[1] = { 512 };

#ifdef __CUDACC__
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
#else
    //assert( (size_t)(allmomenta) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) % mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
#endif

    // Start sigmaKin_lines
    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    allNumerators[ievt] = 0;
    allDenominators[ievt] = 0;
#endif
#else
    const int npagV = nevt / neppV;
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        const unsigned int ievt = ipagV * neppV + ieppV;
        allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        allNumerators[ievt] = 0;
        allDenominators[ievt] = 0;
#endif
      }
    }
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef __CUDACC__
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs );
#endif
#else
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId, nevt );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, nevt );
#endif
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] /= denominators[0]; // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
#endif
#else
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        const unsigned int ievt = ipagV * neppV + ieppV;
        allMEs[ievt] /= denominators[0];                                                 // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
#endif
        //printf( "sigmaKin: ievt=%2d me=%f\n", ievt, allMEs[ievt] );
      }
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

} // end namespace

//==========================================================================
