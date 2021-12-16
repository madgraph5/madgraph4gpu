#ifndef P1_SIGMA_SM_EPEM_MUPMUM_RANDOMNUMBERKERNEL_H 
#define P1_SIGMA_SM_EPEM_MUPMUM_RANDOMNUMBERKERNEL_H 1

#include <future>
#include <vector>

#include "mgOnGpuConfig.h"
//#include "mgOnGpuTypes.h"

// NB This must come AFTER mgOnGpuConfig.h which contains our definition of __global__ when __CUDACC__ is not defined
#ifndef MGONGPU_COMMONRAND_ONHOST
#include "curand.h"
#endif

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  // A base class encapsulating random number generation on a CPU host or on a GPU device
  class RandomNumberKernelBase
  {

  public:

    // === RANDOM NUMBERS ===

    // Get the host buffer[nevt*nparf*4] for the output random numbers (CPU or GPU)
    // [NB on GPU, this is a nullptr unless random numbers are generated on the host]
    fptype* hstRnarray() const { return m_hstRnarray; }

#ifdef __CUDACC__
    // Get the device buffer[nevt*nparf*4] for the output random numbers (GPU)
    fptype* devRnarray() const { return m_devRnarray; }

    // Copy the output random numbers from the host buffer to the device buffer
    // [NB on GPU, this throws unless random numbers are generated on the host]
    void copyHstRnarrayToDevRnArray();
#endif

  protected:

    // Constructor - allocates the output buffer(s) for the given number of events
    RandomNumberKernelBase( int nevt );

    // Destructor - deallocates the output buffer(s)
    virtual ~RandomNumberKernelBase();

  public:

    // Hardcoded parameters (temporarely set them from mgOnGpu; eventually define them only here?)
    static constexpr int nparf = mgOnGpu::nparf;
    static constexpr int np4 = mgOnGpu::np4;

  private:

    // The number of events
    const int m_nevt;

    // The host buffer[nevt*nparf*4] for the output random numbers (CPU or GPU)
    fptype* m_hstRnarray;

#ifdef __CUDACC__
    // The device buffer[nevt*nparf*4] for the output random numbers (GPU)
    fptype* m_devRnarray;
#endif

  };

  //--------------------------------------------------------------------------

  // A class encapsulating common random number generation on a CPU host
  class CommonRandomKernel : public RandomNumberKernelBase
  {
  public:

    // Constructor - allocates the output buffer(s) for the given number of events
    CommonRandomKernel( int nevt );

    // Destructor - deallocates the output buffer(s)
    virtual ~CommonRandomKernel();

    // Seed the random number generator
    //void seedRnGenerator( const int seed );

  private:

    // The common random promises
    std::vector<std::promise<std::vector<fptype>>> m_promises;

  };

  //--------------------------------------------------------------------------

#ifndef MGONGPU_COMMONRAND_ONHOST

  // A class encapsulating CURAND random number generation on a CPU host
  class CurandRandomKernel : public RandomNumberKernelBase
  {
  public:

    // Supported random number generation modes
    enum class RandomNumberMode{ CurandHost=1, CurandDevice=2 };

    // Constructor - allocates the output buffer(s) for the given number of events
    CurandRandomKernel( int nevt );

    // Destructor - deallocates the output buffer(s)
    virtual ~CurandRandomKernel();

    // Seed the random number generator
    void seedRnGenerator( const int seed );

  private:

    // The curand generator
    curandGenerator_t m_rnGen;

  };

#endif

  //--------------------------------------------------------------------------

}
#endif // P1_SIGMA_SM_EPEM_MUPMUM_RANDOMNUMBERKERNEL_H
