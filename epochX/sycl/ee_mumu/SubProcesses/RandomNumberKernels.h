#ifndef RANDOMNUMBERKERNELS_H 
#define RANDOMNUMBERKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "MemoryBuffers.h"

namespace mg5amcGpu
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
    virtual void seedGenerator( const int seed ) = 0;

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
    RandomNumberKernelBase( BufferRandomNumbers& rnarray ) : m_rnarray( rnarray ){}

  public:

    // Destructor
    virtual ~RandomNumberKernelBase(){}
    
    // Seed the random number generator
    virtual void seedGenerator( const int seed ) = 0;

    // Generate the random number array
    virtual void generateRnarray() = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The buffer for the output random numbers
    BufferRandomNumbers& m_rnarray;

  };

  //--------------------------------------------------------------------------

  // A class encapsulating common random number generation on a CPU host
  class CommonRandomNumberKernel final : public RandomNumberKernelBase
  {
  public:

    // Constructor from an existing output buffer
    CommonRandomNumberKernel( BufferRandomNumbers& rnarray );

    // Destructor
    ~CommonRandomNumberKernel(){}

    // Seed the random number generator
    void seedGenerator( const int seed ) override final { m_seed = seed; };

    // Generate the random number array
    void generateRnarray() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  private:

    // The generator seed
    int m_seed;

  };

  //--------------------------------------------------------------------------

}
#endif // RANDOMNUMBERKERNELS_H
