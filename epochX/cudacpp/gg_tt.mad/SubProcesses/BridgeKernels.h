#ifndef BRIDGEKERNELS_H
#define BRIDGEKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "Bridge.h"
#include "MatrixElementKernels.h"
#include "MemoryBuffers.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  // A Bridge wrapper base class encapsulating matrix element calculations on a CPU host
  class BridgeKernelBase : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    BridgeKernelBase( const BufferMomenta& momenta,         // input: momenta
                      BufferMatrixElements& matrixElements, // output: matrix elements
                      const size_t nevt );

    // Destructor
    virtual ~BridgeKernelBase() {}

    // Transpose input momenta from C to Fortran before the matrix element calculation in the Bridge
    virtual void transposeInputMomentaC2F() = 0;

  protected:

    // The wrapped bridge
    Bridge<fptype> m_bridge;
  };

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A Bridge wrapper class encapsulating matrix element calculations on a CPU host
  class BridgeKernelHost final : public BridgeKernelBase
  {
  public:

    // Constructor from existing input and output buffers
    BridgeKernelHost( const BufferMomenta& momenta,         // input: momenta
                      BufferMatrixElements& matrixElements, // output: matrix elements
                      const size_t nevt );

    // Destructor
    virtual ~BridgeKernelHost() {}

    // Transpose input momenta from C to Fortran before the matrix element calculation in the Bridge
    void transposeInputMomentaC2F() override final;

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  private:

    // The buffer for the input momenta, transposed to Fortran array indexing
    HostBufferMomenta m_fortranMomenta;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A Bridge wrapper class encapsulating matrix element calculations on a GPU device
  class BridgeKernelDevice : public BridgeKernelBase
  {
  public:

    // Constructor from existing input and output buffers
    BridgeKernelDevice( const BufferMomenta& momenta,         // input: momenta
                        BufferMatrixElements& matrixElements, // output: matrix elements
                        const size_t gpublocks,
                        const size_t gputhreads );

    // Destructor
    virtual ~BridgeKernelDevice() {}

    // Transpose input momenta from C to Fortran before the matrix element calculation in the Bridge
    void transposeInputMomentaC2F() override final;

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The buffer for the input momenta, transposed to Fortran array indexing
    PinnedHostBufferMomenta m_fortranMomenta;

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;
  };
#endif

  //--------------------------------------------------------------------------
}
#endif // BRIDGEKERNELS_H
