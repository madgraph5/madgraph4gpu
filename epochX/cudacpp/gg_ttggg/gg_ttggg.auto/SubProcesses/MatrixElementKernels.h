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
    MatrixElementKernelBase( const BufferMomenta& momenta,          // input: momenta
                             BufferMatrixElements& matrixElements ) // output: matrix elements
      : m_momenta( momenta )
      , m_matrixElements( matrixElements ){}

  public:

    // Destructor
    virtual ~MatrixElementKernelBase(){}

    // Compute good helicities
    virtual void computeGoodHelicities() = 0;

    // Compute matrix elements
    virtual void computeMatrixElements() = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The buffer for the input momenta
    const BufferMomenta& m_momenta;

    // The buffer for the output matrix elements
    BufferMatrixElements& m_matrixElements;

  };

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A class encapsulating matrix element calculations on a CPU host
  class MatrixElementKernelHost final : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelHost( const BufferMomenta& momenta,         // input: momenta
                             BufferMatrixElements& matrixElements, // output: matrix elements
                             const size_t nevt );

    // Destructor
    virtual ~MatrixElementKernelHost(){}

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

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
                               BufferMatrixElements& matrixElements, // output: matrix elements
                               const size_t gpublocks,
                               const size_t gputhreads );

    // Destructor
    virtual ~MatrixElementKernelDevice(){}

    // Reset gpublocks and gputhreads
    void setGrid( const int gpublocks, const int gputhreads );

    // Compute good helicities
    void computeGoodHelicities() override final;

    // Compute matrix elements
    void computeMatrixElements() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;

  };
#endif

  //--------------------------------------------------------------------------

}
#endif // MATRIXELEMENTKERNELS_H
