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
// Process: g g > t t~ g WEIGHTED<=3 @1

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
    constexpr int ncolor = 6;

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

      // *** DIAGRAM 1 OF 16 ***

      // Wavefunction(s) for diagram number 1
      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][0], -1, w_fp[0], 0 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][1], -1, w_fp[1], 1 );

      oxxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][2], +1, w_fp[2], 2 );

      ixxxxx<M_ACCESS, W_ACCESS>( momenta, cIPD[0], cHel[ihel][3], -1, w_fp[3], 3 );

      vxxxxx<M_ACCESS, W_ACCESS>( momenta, 0., cHel[ihel][4], +1, w_fp[4], 4 );

      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], COUPs[0], 0., 0., w_fp[5] );
      FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], COUPs[1], 0., 0., w_fp[6] );

      // Amplitude(s) for diagram number 1
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[6], w_fp[4], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 1 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] -= amp_sv[0];
      jamp_sv[2] += amp_sv[0];
      jamp_sv[4] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];

      // *** DIAGRAM 2 OF 16 ***

      // Wavefunction(s) for diagram number 2
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[7] );

      // Amplitude(s) for diagram number 2
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[7], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 2 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 3 OF 16 ***

      // Wavefunction(s) for diagram number 3
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[4], COUPs[1], cIPD[0], cIPD[1], w_fp[8] );

      // Amplitude(s) for diagram number 3
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[2], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 3 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 4 OF 16 ***

      // Wavefunction(s) for diagram number 4
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[5] );
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[9] );

      // Amplitude(s) for diagram number 4
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[5], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 4 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] -= amp_sv[0];

      // *** DIAGRAM 5 OF 16 ***

      // Wavefunction(s) for diagram number 5
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[0], 0., 0., w_fp[10] );

      // Amplitude(s) for diagram number 5
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[5], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 5 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 6 OF 16 ***

      // Wavefunction(s) for diagram number 6
      // (none)

      // Amplitude(s) for diagram number 6
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 6 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] -= amp_sv[0];

      // *** DIAGRAM 7 OF 16 ***

      // Wavefunction(s) for diagram number 7
      FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], cIPD[0], cIPD[1], w_fp[5] );
      FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], cIPD[0], cIPD[1], w_fp[11] );

      // Amplitude(s) for diagram number 7
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[11], w_fp[4], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 7 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] -= amp_sv[0];

      // *** DIAGRAM 8 OF 16 ***

      // Wavefunction(s) for diagram number 8
      // (none)

      // Amplitude(s) for diagram number 8
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[2], w_fp[10], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 8 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 9 OF 16 ***

      // Wavefunction(s) for diagram number 9
      // (none)

      // Amplitude(s) for diagram number 9
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[7], w_fp[1], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 9 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[5] -= amp_sv[0];

      // *** DIAGRAM 10 OF 16 ***

      // Wavefunction(s) for diagram number 10
      VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], COUPs[0], 0., 0., w_fp[5] );

      // Amplitude(s) for diagram number 10
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[11], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 10 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 11 OF 16 ***

      // Wavefunction(s) for diagram number 11
      // (none)

      // Amplitude(s) for diagram number 11
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[2], w_fp[5], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 11 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
      jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];

      // *** DIAGRAM 12 OF 16 ***

      // Wavefunction(s) for diagram number 12
      // (none)

      // Amplitude(s) for diagram number 12
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[5], w_fp[1], w_fp[6], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 12 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[1] += amp_sv[0];
      jamp_sv[2] -= amp_sv[0];
      jamp_sv[3] += amp_sv[0];
      jamp_sv[4] -= amp_sv[0];

      // *** DIAGRAM 13 OF 16 ***

      // Wavefunction(s) for diagram number 13
      // (none)

      // Amplitude(s) for diagram number 13
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[11], w_fp[0], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 13 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[2] -= amp_sv[0];

      // *** DIAGRAM 14 OF 16 ***

      // Wavefunction(s) for diagram number 14
      // (none)

      // Amplitude(s) for diagram number 14
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[9], w_fp[7], w_fp[0], COUPs[1], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 14 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[4] -= amp_sv[0];

      // *** DIAGRAM 15 OF 16 ***

      // Wavefunction(s) for diagram number 15
      // (none)

      // Amplitude(s) for diagram number 15
      VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[10], w_fp[6], COUPs[0], &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId == 15 ) numerators_sv += cxabs2( amp_sv[0] );
      if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
      jamp_sv[0] += amp_sv[0];
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[5] += amp_sv[0];

      // *** DIAGRAM 16 OF 16 ***

      // Wavefunction(s) for diagram number 16
      VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[10] );
      VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[6] );
      VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 0., 0., w_fp[9] );

      // Amplitude(s) for diagram number 16
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], w_fp[10], COUPs[1], &amp_fp[0] );
      jamp_sv[0] += amp_sv[0];
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[5] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], w_fp[6], COUPs[1], &amp_fp[0] );
      jamp_sv[1] -= amp_sv[0];
      jamp_sv[2] += amp_sv[0];
      jamp_sv[3] -= amp_sv[0];
      jamp_sv[4] += amp_sv[0];
      FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[2], w_fp[9], COUPs[1], &amp_fp[0] );
      jamp_sv[0] -= amp_sv[0];
      jamp_sv[2] += amp_sv[0];
      jamp_sv[4] += amp_sv[0];
      jamp_sv[5] -= amp_sv[0];

      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called CPPProcess::matrix_1_gg_ttxg()?)

      // The color denominators (initialize all array elements, with ncolor=6)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype denom[ncolor] = { 9, 9, 9, 9, 9, 9 }; // 1-D array[6]

      // The color matrix (initialize all array elements, with ncolor=6)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype cf[ncolor][ncolor] = {
        { 64, -8, -8, 1, 1, 10 },
        { -8, 64, 1, 10, -8, 1 },
        { -8, 1, 64, -8, 10, 1 },
        { 1, 10, -8, 64, 1, -8 },
        { 1, -8, 10, 1, 64, -8 },
        { 10, 1, 1, -8, -8, 64 } }; // 2-D array[6][6]

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
      for( int icol = 0; icol < ncolor; icol++ )
      {
        cxtype_sv ztemp_sv = cxzero_sv();
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp_sv += cf[icol][jcol] * jamp_sv[jcol];
        // OLD implementation: why is this not slower? maybe compiler does not compute imaginary part of "ztemp_sv*cxconj(jamp_sv[icol])"?
        //deltaMEs += cxreal( ztemp_sv * cxconj( jamp_sv[icol] ) ) / denom[icol];
        // NEW implementation: keep this even if (surprisingly) it is not faster! it is clearer and easier for tensor core offload anyway...
        // Rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB!
        deltaMEs += ( cxreal( ztemp_sv ) * cxreal( jamp_sv[icol] ) + cximag( ztemp_sv ) * cximag( jamp_sv[icol] ) ) / denom[icol];
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
      { -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, 1 },
      { -1, -1, -1, 1, -1 },
      { -1, -1, -1, 1, 1 },
      { -1, -1, 1, -1, -1 },
      { -1, -1, 1, -1, 1 },
      { -1, -1, 1, 1, -1 },
      { -1, -1, 1, 1, 1 },
      { -1, 1, -1, -1, -1 },
      { -1, 1, -1, -1, 1 },
      { -1, 1, -1, 1, -1 },
      { -1, 1, -1, 1, 1 },
      { -1, 1, 1, -1, -1 },
      { -1, 1, 1, -1, 1 },
      { -1, 1, 1, 1, -1 },
      { -1, 1, 1, 1, 1 },
      { 1, -1, -1, -1, -1 },
      { 1, -1, -1, -1, 1 },
      { 1, -1, -1, 1, -1 },
      { 1, -1, -1, 1, 1 },
      { 1, -1, 1, -1, -1 },
      { 1, -1, 1, -1, 1 },
      { 1, -1, 1, 1, -1 },
      { 1, -1, 1, 1, 1 },
      { 1, 1, -1, -1, -1 },
      { 1, 1, -1, -1, 1 },
      { 1, 1, -1, 1, -1 },
      { 1, 1, -1, 1, 1 },
      { 1, 1, 1, -1, -1 },
      { 1, 1, 1, -1, 1 },
      { 1, 1, 1, 1, -1 },
      { 1, 1, 1, 1, 1 } };
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
    constexpr int denominators[1] = { 256 };

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
