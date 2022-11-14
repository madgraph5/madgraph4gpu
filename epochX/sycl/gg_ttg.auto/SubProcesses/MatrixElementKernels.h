#ifndef MATRIXELEMENTKERNELS_H 
#define MATRIXELEMENTKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "MemoryBuffers.h"

namespace mg5amcGpu
{

  //--------------------------------------------------------------------------

  // A base class encapsulating matrix element calculations on a CPU host or on a GPU device
  class MatrixElementKernelBase //: virtual public IMatrixElementKernel
  {
  protected:

    // Constructor from existing input and output buffers
    MatrixElementKernelBase( const BufferMomenta& momenta,          // input: momenta
                             BufferMatrixElements& matrixElements,  // output: matrix elements
                             sycl::queue q )
      : m_momenta( momenta )
      , m_matrixElements( matrixElements )
      , m_q( q ){}

  public:

    // Destructor
    virtual ~MatrixElementKernelBase(){}

    // Set physics parameters
    virtual void setDeviceArrays( const short* tHel, const cxtype* tIPC, const fptype* tIPD ) = 0;

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
    
    // The sycl::queue to run kernels on
    sycl::queue m_q;

  };

  //--------------------------------------------------------------------------

  // A class encapsulating matrix element calculations on a GPU device
  class MatrixElementKernelDevice : public MatrixElementKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    MatrixElementKernelDevice( const BufferMomenta& momenta,         // input: momenta
                               BufferMatrixElements& matrixElements, // output: matrix elements
                               sycl::queue q,
                               const size_t gpublocks,
                               const size_t gputhreads );

    // Destructor
    virtual ~MatrixElementKernelDevice(){}

    // Reset gpublocks and gputhreads
    void setGrid( const int gpublocks, const int gputhreads );

    // Set physics parameters
    void setDeviceArrays( const short* tHel, const cxtype* tIPC, const fptype* tIPD ) override final;

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

    // FIXME add comment about variables below (maybe typedef in MemoryBuffer instead of DeviceBufferBase)
    DeviceBufferBase<short> m_cHel;
    DeviceBufferBase<fptype> m_cIPC;
    DeviceBufferBase<fptype> m_cIPD;

    DeviceBufferBase<int> m_cNGoodHel;
    DeviceBufferBase<int> m_cGoodHel;
  };

  //--------------------------------------------------------------------------

}
#endif // MATRIXELEMENTKERNELS_H
