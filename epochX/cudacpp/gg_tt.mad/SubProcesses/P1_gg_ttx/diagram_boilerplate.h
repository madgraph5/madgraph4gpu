// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Nov 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024-2025) for the MG5aMC CUDACPP plugin.

#ifdef MGONGPUCPP_GPUIMPL

    //-------------
    // GPU only
    //-------------

    //using namespace mg5amcGpu;
    using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using CI_ACCESS = DeviceAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
    using J_ACCESS = DeviceAccessJamp;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
    // SCALAR channelId for the current event (CUDA)
    unsigned int channelId = gpu_channelId( channelIds );
#endif

    // Wavefunctions
    // Buffer wfs for one helicity and nevt events is a DeviceBufferSimple with ( nwf * nevt * nw6 * nx2 ) fptypes
    // The striding between the nwf wavefunction buffers is ( nevt * nw6 * nx2 ) fptypes
    // Internally diagramXXX methods pass a w_fp[iwf] to ixx/FFV methods (as argument 'fptype wavefunctions[]')
    // Internally ixx/FFV methods call 'cxtype_sv* fi = W_ACCESS::kernelAccess( wavefunctions )' and then use fi[iw6]
    // This means that the fi pointer must point to a [RIRIRIRIRIRI] contiguous buffer of size nw6*nx2=12
    // The striding between events is nw6*nx2=12 and this is what W_ACCESS::kernelAccess must respect
    // (En passant, note that this means that events cannot be contiguous in the present code, memory is not coalesced)
    const int nevt = gridDim.x * blockDim.x;
    fptype* w_fp[nwf];
    for( int iwf = 0; iwf < nwf; iwf++ ) w_fp[iwf] = wfs + iwf * nevt * nw6 * mgOnGpu::nx2;

    // Couplings
    constexpr size_t nxcoup = ndcoup + nIPC; // both dependent and independent couplings (FIX #823: nIPC instead of nicoup)
    const fptype* allCOUPs[nxcoup];
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 186 // e.g. <<warning #186-D: pointless comparison of unsigned integer with zero>>
#endif
    // Dependent couplings, vary event-by-event
    for( size_t idcoup = 0; idcoup < ndcoup; idcoup++ )
      allCOUPs[idcoup] = CD_ACCESS::idcoupAccessBufferConst( couplings, idcoup );
    // Independent couplings, fixed for all events
    for( size_t iicoup = 0; iicoup < nIPC; iicoup++ ) // (FIX #823: nIPC instead of nicoup)
      allCOUPs[ndcoup + iicoup] = CI_ACCESS::iicoupAccessBufferConst( cIPC, iicoup );
#ifdef __CUDACC__ // this must be __CUDACC__ (not MGONGPUCPP_GPUIMPL)
#pragma nv_diagnostic pop
#endif
    const fptype* COUPs[nxcoup];
    for( size_t ixcoup = 0; ixcoup < nxcoup; ixcoup++ ) COUPs[ixcoup] = allCOUPs[ixcoup];

#else

    //-------------
    // C++ only
    //-------------

    //using namespace mg5amcCpu;
    using W_ACCESS = HostAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = HostAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using J_ACCESS = HostAccessJamp;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
    // SCALAR channelId for the current SIMD event page (C++)
    unsigned int channelId = *channelIds;
#endif

    // Wavefunctions
    // Reinterpret wfs as "cxtype_sv w_sv[nwf][nw6]" and build "fptype* w_fp[nwf]" where "w_fp[iwf] = (fptype*)( w_sv[iwf] )"
    fptype (*w_fp)[nw6*sizeof(cxtype_sv)/sizeof(fptype)] = (fptype (*)[nw6*sizeof(cxtype_sv)/sizeof(fptype)])(wfs);

#endif

    //-------------
    // GPU or C++
    //-------------

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Numerators and denominators for the current event (CUDA) or SIMD event page (C++)
    fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
    fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
#else
    assert( channelIds == nullptr );
    assert( numerators == nullptr );
    assert( denominators == nullptr );
#endif

    // Local TEMPORARY variables for a subset of Feynman diagrams in the given CUDA event (ievt) or C++ event page (ipagV)
    cxtype_sv amp_sv[1]; // invariant amplitude for one given Feynman diagram
    fptype* amp_fp;      // proof of concept for using fptype* in the interface
    amp_fp = reinterpret_cast<fptype*>( amp_sv );
