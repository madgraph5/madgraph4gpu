// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Nov 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#ifdef MGONGPUCPP_GPUIMPL

    //-------------
    // GPU only
    //-------------

    //using namespace mg5amcGpu;
    using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
    using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
    using J_ACCESS = DeviceAccessJamp;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    using NUM_ACCESS = DeviceAccessNumerators;    // non-trivial access: buffer includes all events
    using DEN_ACCESS = DeviceAccessDenominators;  // non-trivial access: buffer includes all events
    // SCALAR channelId for the current event (CUDA)
    unsigned int channelId = gpu_channelId( channelIds );
#endif

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
