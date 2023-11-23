// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.

#ifndef MATRIXELEMENTKERNELS_H
#define MATRIXELEMENTKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "MemoryBuffers.h"

#ifdef __CUDACC__
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
                             BufferMatrixElements& matrixElements, // output: matrix elements
                             BufferSelectedHelicity& selhel,       // output: helicity selection
                             BufferSelectedColor& selcol )         // output: color selection
      : m_momenta( momenta )
      , m_gs( gs )
      , m_rndhel( rndhel )
      , m_rndcol( rndcol )
      , m_matrixElements( matrixElements )
      , m_selhel( selhel )
      , m_selcol( selcol )
    {
    }

  public:

    // Destructor
    virtual ~MatrixElementKernelBase() {}

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    virtual int computeGoodHelicities() = 0;

    // Compute matrix elements
    virtual void computeMatrixElements( const unsigned int channelId ) = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The buffer for the input momenta
    const BufferMomenta& m_momenta;

    // The buffer for the gs to calculate the alphaS values
    const BufferGs& m_gs;

    // The buffer for the random numbers for helicity selection
    const BufferRndNumHelicity& m_rndhel;

    // The buffer for the random numbers for color selection
    const BufferRndNumColor& m_rndcol;

    // The buffer for the output matrix elements
    BufferMatrixElements& m_matrixElements;

    // The buffer for the output helicity selection
    BufferSelectedHelicity& m_selhel;

    // The buffer for the output color selection
    BufferSelectedColor& m_selcol;
  };

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A class encapsulating matrix element calculations on a CPU host
  class MatrixElementKernelHost final : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                             const BufferGs& gs,                   // input: gs for alphaS
                             const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                             const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                             BufferMatrixElements& matrixElements, // output: matrix elements
                             BufferSelectedHelicity& selhel,       // output: helicity selection
                             BufferSelectedColor& selcol,          // output: color selection
                             const size_t nevt );

    // Destructor
    virtual ~MatrixElementKernelHost() {}

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    int computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements( const unsigned int channelId ) override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

    // Does this host system support the SIMD used in the matrix element calculation?
    // [NB: SIMD vectorization in mg5amc C++ code is currently only used in the ME calculations below MatrixElementKernelHost!]
    static bool hostSupportsSIMD( const bool verbose = true );

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

#ifdef __CUDACC__
  // A class encapsulating matrix element calculations on a GPU device
  class MatrixElementKernelDevice : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                               const BufferGs& gs,                   // input: gs for alphaS
                               const BufferRndNumHelicity& rndhel,   // input: random numbers for helicity selection
                               const BufferRndNumColor& rndcol,      // input: random numbers for color selection
                               BufferMatrixElements& matrixElements, // output: matrix elements
                               BufferSelectedHelicity& selhel,       // output: helicity selection
                               BufferSelectedColor& selcol,          // output: color selection
                               const size_t gpublocks,
                               const size_t gputhreads );

    // Destructor
    virtual ~MatrixElementKernelDevice() {}

    // Reset gpublocks and gputhreads
    void setGrid( const int gpublocks, const int gputhreads );

    // Compute good helicities (returns nGoodHel, the number of good helicity combinations out of ncomb)
    int computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements( const unsigned int channelId ) override final;

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

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;
  };
#endif

  //--------------------------------------------------------------------------
}
#endif // MATRIXELEMENTKERNELS_H
