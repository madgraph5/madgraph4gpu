#ifndef RANDOMNUMBERKERNELS_H
#define RANDOMNUMBERKERNELS_H 1

#include "mgOnGpuConfig.h"

// NB This must come AFTER mgOnGpuConfig.h which contains our definition of __global__ when __CUDACC__ is not defined
#ifndef MGONGPU_HAS_NO_CURAND
#include "curand.h"
#endif

#include "MemoryBuffers.h"

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  /*
  // An interface encapsulating random number generation on a CPU host or on a GPU device
  class IRandomNumberKernel
  {
  public:

    // Destructor
    virtual ~IRandomNumberKernel(){}

    // Seed the random number generator
    virtual void seedGenerator( const unsigned int seed ) = 0;

    // Generate the random number array
    virtual void generateRnarray() = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  };
  */

  //--------------------------------------------------------------------------

  // A base class encapsulating random number generation on a CPU host or on a GPU device
  class RandomNumberKernelBase //: virtual public IRandomNumberKernel
  {

  protected:

    // Constructor from an existing output buffer
    RandomNumberKernelBase( BufferRndNumMomenta& rnarray )
      : m_rnarray( rnarray ) {}

  public:

    // Destructor
    virtual ~RandomNumberKernelBase() {}

    // Seed the random number generator
    virtual void seedGenerator( const unsigned int seed ) = 0;

    // Generate the random number array
    virtual void generateRnarray() = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The buffer for the output random numbers
    BufferRndNumMomenta& m_rnarray;
  };

  //--------------------------------------------------------------------------

  // A class encapsulating common random number generation on a CPU host
  class CommonRandomNumberKernel final : public RandomNumberKernelBase
  {
  public:

    // Constructor from an existing output buffer
    CommonRandomNumberKernel( BufferRndNumMomenta& rnarray );

    // Destructor
    ~CommonRandomNumberKernel() {}

    // Seed the random number generator
    void seedGenerator( const unsigned int seed ) override final { m_seed = seed; };

    // Generate the random number array
    void generateRnarray() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  private:

    // The generator seed
    unsigned int m_seed;
  };

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HAS_NO_CURAND
  // A class encapsulating CURAND random number generation on a CPU host or on a GPU device
  class CurandRandomNumberKernel final : public RandomNumberKernelBase
  {
  public:

    // Constructor from an existing output buffer
    CurandRandomNumberKernel( BufferRndNumMomenta& rnarray, const bool onDevice );

    // Destructor
    ~CurandRandomNumberKernel();

    // Seed the random number generator
    void seedGenerator( const unsigned int seed ) override final;

    // Generate the random number array
    void generateRnarray() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return m_isOnDevice; }

  private:

    // Create the generator (workaround for #429: do this in every seedGenerator call rather than only in the ctor)
    void createGenerator();

    // Destroy the generator (workaround for #429: do this in every seedGenerator call rather than only in the ctor)
    void destroyGenerator();

  private:

    // Is this a host or device kernel?
    const bool m_isOnDevice;

    // The curand generator
    curandGenerator_t m_rnGen;
  };

#endif

  //--------------------------------------------------------------------------
}
#endif // RANDOMNUMBERKERNELS_H
