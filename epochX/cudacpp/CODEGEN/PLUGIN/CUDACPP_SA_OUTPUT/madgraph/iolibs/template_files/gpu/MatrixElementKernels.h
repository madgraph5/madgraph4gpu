// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef MATRIXELEMENTKERNELS_H
#define MATRIXELEMENTKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "MemoryBuffers.h"

#include <map>

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  // A base class encapsulating matrix element calculations on a CPU host or on a GPU device
  class MatrixElementKernelBase //: virtual public IMatrixElementKernel
  {
  protected:

    // Constructor from existing input and output buffers
    MatrixElementKernelBase( const BufferMomenta& momenta,         // input: momenta
                             const BufferGs& gs,                   // input: gs for alphaS
                             const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                             const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                             const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                             BufferMatrixElements& matrixElements, // output: matrix elements
                             BufferSelectedHelicity& selhel,       // output: helicity selection
                             BufferSelectedColor& selcol );        // output: color selection

  public:

    // Destructor
    virtual ~MatrixElementKernelBase();

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    virtual int computeGoodHelicities() = 0;

    // Compute matrix elements
    virtual void computeMatrixElements( const bool useChannelIds ) = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

    // Dump signalling FPEs (#831 and #837)
    static void dumpSignallingFPEs();

#ifdef MGONGPU_CHANNELID_DEBUG
    // Add a MEK identifier for the channelId debug printout
    void setTagForNevtProcessedByChannel( const std::string& tag ) { m_tag = tag; }

  protected:
    // Update number of events processed by channel
    void updateNevtProcessedByChannel( const unsigned int* pHstChannelIds, const size_t nevt );

    // Dump number of events processed by channel
    void dumpNevtProcessedByChannel();
#endif

  protected:

    // The buffer for the input momenta
    const BufferMomenta& m_momenta;

    // The buffer for the gs to calculate the alphaS values
    const BufferGs& m_gs;

    // The buffer for the random numbers for helicity selection
    const BufferRndNumHelicity& m_rndhel;

    // The buffer for the random numbers for color selection
    const BufferRndNumColor& m_rndcol;

    // The buffer for the channel ids for single-diagram enhancement
    const BufferChannelIds& m_channelIds;

    // The buffer for the output matrix elements
    BufferMatrixElements& m_matrixElements;

    // The buffer for the output helicity selection
    BufferSelectedHelicity& m_selhel;

    // The buffer for the output color selection
    BufferSelectedColor& m_selcol;

#ifdef MGONGPU_CHANNELID_DEBUG
    // The events-per-channel counter for debugging
    std::map<size_t, size_t> m_nevtProcessedByChannel;

    // The tag for events-per-channel debugging
    std::string m_tag;
#endif
  };

  //--------------------------------------------------------------------------

#ifndef MGONGPUCPP_GPUIMPL
  // A class encapsulating matrix element calculations on a CPU host
  class MatrixElementKernelHost final : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                             const BufferGs& gs,                   // input: gs for alphaS
                             const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                             const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                             const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                             BufferMatrixElements& matrixElements, // output: matrix elements
                             BufferSelectedHelicity& selhel,       // output: helicity selection
                             BufferSelectedColor& selcol,          // output: color selection
                             const size_t nevt );

    // Destructor
    virtual ~MatrixElementKernelHost();

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    int computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements( const bool useChannelIds ) override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  private:

    // Does this host system support the SIMD used in the matrix element calculation?
    // [NB: this is private, SIMD vectorization in mg5amc C++ code is currently only used in the ME calculations below MatrixElementKernelHost!]
    static bool hostSupportsSIMD( const bool verbose = false ); // ZW: default verbose false

  private:

    // The buffer for the event-by-event couplings that depends on alphas QCD
    HostBufferCouplings m_couplings;

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // The buffer for the event-by-event numerators of multichannel factors
    HostBufferNumerators m_numerators;

    // The buffer for the event-by-event denominators of multichannel factors
    HostBufferDenominators m_denominators;
#endif
  };
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  // A class encapsulating matrix element calculations on a GPU device
  class MatrixElementKernelDevice : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                               const BufferGs& gs,                   // input: gs for alphaS
                               const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                               const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                               const BufferChannelIds& channelIds,   // input: channel ids for single-diagram enhancement
                               BufferMatrixElements& matrixElements, // output: matrix elements
                               BufferSelectedHelicity& selhel,       // output: helicity selection
                               BufferSelectedColor& selcol,          // output: color selection
                               const size_t gpublocks,
                               const size_t gputhreads );

    // Destructor
    virtual ~MatrixElementKernelDevice();

    // Reset gpublocks and gputhreads
    void setGrid( const int gpublocks, const int gputhreads );

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    int computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements( const bool useChannelIds ) override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The buffer for the event-by-event couplings that depends on alphas QCD
    DeviceBufferCouplings m_couplings;

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // The buffer for the event-by-event numerators of multichannel factors
    DeviceBufferNumerators m_numerators;

    // The buffer for the event-by-event denominators of multichannel factors
    DeviceBufferDenominators m_denominators;
#endif

#ifdef MGONGPU_CHANNELID_DEBUG
    // The **host** buffer for the channelId array
    // FIXME? MEKD should accept a host buffer as an argument instead of a device buffer, so that a second copy can be avoided?
    PinnedHostBufferChannelIds m_hstChannelIds;
#endif

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;
  };
#endif

  //--------------------------------------------------------------------------
}
#endif // MATRIXELEMENTKERNELS_H
