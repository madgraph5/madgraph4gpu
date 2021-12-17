#ifndef RANDOMNUMBERKERNELS_H 
#define RANDOMNUMBERKERNELS_H 1

#include "mgOnGpuConfig.h"

// NB This must come AFTER mgOnGpuConfig.h which contains our definition of __global__ when __CUDACC__ is not defined
#ifndef MGONGPU_HAS_NO_CURAND
#include "curand.h"
#endif

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  // Supported random number generation modes
  enum class RandomNumberMode{ CommonRandom=0, CurandHost=1, CurandDevice=2 };

  //--------------------------------------------------------------------------

  // A base class encapsulating random number generation on a CPU host or on a GPU device
  class RandomNumberKernelBase
  {

  protected:

#ifndef __CUDACC__
    // Constructor - allocates the output buffer(s) for the given number of events
    RandomNumberKernelBase( const int nevt );
#else
    // Constructor - allocates the output buffer(s) for the given number of events
    RandomNumberKernelBase( const int nevt, const bool useHstRnarray = true );
#endif

  public:

    // Destructor - deallocates the output buffer(s)
    virtual ~RandomNumberKernelBase();

    // Seed the random number generator (throws if seed is <= 0)
    virtual void seedGenerator( const int seed );

    // Generate the random number array (throws if seed is <= 0)
    virtual void generateRnarray();

    // === RANDOM NUMBERS ===

    // Get the host buffer[nevt*nparf*4] for the output random numbers (CPU or GPU)
    // [NB on GPU, this is a nullptr unless random numbers are generated on the host]
    fptype* hstRnarray() const { return m_hstRnarray; }

#ifdef __CUDACC__
    // Get the device buffer[nevt*nparf*4] for the output random numbers (GPU)
    fptype* devRnarray() const { return m_devRnarray; }

    // Copy the output random numbers from the host buffer to the device buffer
    // [NB on GPU, this throws unless random numbers are generated on the host]
    void copyHstRnarrayToDevRnarray();
#endif

    // Get the number of elements in the random number buffer(s)
    int nRnarray() const { return np4 * nparf * m_nevt; }

  public:

    // Hardcoded parameters (temporarely set them from mgOnGpu; eventually define them only here?)
    static constexpr int nparf = mgOnGpu::nparf;
    static constexpr int np4 = mgOnGpu::np4;
#ifndef __CUDACC__
    static constexpr int cppAlign = mgOnGpu::cppAlign;
#endif

  protected:

    // The number of events
    const int m_nevt;

    // The generator seed
    int m_seed;

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
    CommonRandomKernel( const int nevt );

    // Destructor - deallocates the output buffer(s)
    virtual ~CommonRandomKernel();

    // Generate the random number array (throws if seed is <= 0)
    void generateRnarray();

  };

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HAS_NO_CURAND

  // A class encapsulating CURAND random number generation on a CPU host
  class CurandRandomKernel : public RandomNumberKernelBase
  {
  public:

    // Constructor - allocates the output buffer(s) for the given number of events
    CurandRandomKernel( int nevt, RandomNumberMode mode );

    // Destructor - deallocates the output buffer(s)
    virtual ~CurandRandomKernel();

    // Seed the random number generator
    void seedGenerator( const int seed );

    // Generate the random number array (throws if seed is <= 0)
    void generateRnarray();

  private:

    // The random number generation mode
    const RandomNumberMode m_mode;

    // The curand generator
    curandGenerator_t m_rnGen;

  };

#endif

  //--------------------------------------------------------------------------

}
#endif // RANDOMNUMBERKERNELS_H
