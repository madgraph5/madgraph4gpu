#ifndef BRIDGEKERNELS_H 
#define BRIDGEKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "Bridge.h"
#include "MatrixElementKernels.h"
#include "MemoryBuffers.h"

#undef BRIDGEDEBUG
#define BRIDGEDEBUG 1 // DEBUG!

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A Bridge wrapper class encapsulating matrix element calculations on a CPU host 
  class BridgeKernelHost final : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    BridgeKernelHost( const BufferMomenta& momenta,         // input: momenta
                      BufferMatrixElements& matrixElements, // output: matrix elements
                      const size_t nevt );
    
    // Destructor
    virtual ~BridgeKernelHost(){}

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  private:

    // The wrapped bridge
    Bridge<fptype> m_bridge;

    // The buffer for the input momenta, transposed to Fortran array indexing
    HostBufferMomenta m_fortranMomenta;

#ifdef BRIDGEDEBUG
    HostBufferMomenta m_momenta2;
    MatrixElementKernelHost m_mek;
#endif

  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A Bridge wrapper class encapsulating matrix element calculations on a GPU device
  class BridgeKernelDevice : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    BridgeKernelDevice( const BufferMomenta& momenta,         // input: momenta
                        BufferMatrixElements& matrixElements, // output: matrix elements
                        const size_t gpublocks,
                        const size_t gputhreads );

    // Destructor
    virtual ~BridgeKernelDevice(){}

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The wrapped bridge
    Bridge<fptype> m_bridge;
    
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
