#ifndef CROSSSECTIONKERNELS_H 
#define CROSSSECTIONKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "EventStatistics.h"
#include "MemoryBuffers.h"

//============================================================================

// Disabling fast math is essential here, otherwise results are undefined
// See https://stackoverflow.com/a/40702790 about __attribute__ on gcc
// See https://stackoverflow.com/a/32292725 about __attribute__ on clang
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline bool fp_is_abnormal( const fptype& fp )
{
  if ( std::isnan( fp ) ) return true;
  if ( fp != fp ) return true;
  return false;
}

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline bool fp_is_zero( const fptype& fp )
{
  if ( fp == 0 ) return true;
  return false;
}

// See https://en.cppreference.com/w/cpp/numeric/math/FP_categories
#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline const char* fp_show_class( const fptype& fp )
{
  switch( std::fpclassify( fp ) ) {
  case FP_INFINITE:  return "Inf";
  case FP_NAN:       return "NaN";
  case FP_NORMAL:    return "normal";
  case FP_SUBNORMAL: return "subnormal";
  case FP_ZERO:      return "zero";
  default:           return "unknown";
  }
}

#ifdef __clang__
__attribute__((optnone))
#else
__attribute__((optimize("-fno-fast-math")))
#endif
inline void debug_me_is_abnormal( const fptype& me, size_t ievtALL )
{
  std::cout << "DEBUG[" << ievtALL << "]"
            << " ME=" << me
            << " fpisabnormal=" << fp_is_abnormal( me )
            << " fpclass=" << fp_show_class( me )
            << " (me==me)=" << ( me == me )
            << " (me==me+1)=" << ( me == me+1 )
            << " isnan=" << std::isnan( me )
            << " isfinite=" << std::isfinite( me )
            << " isnormal=" << std::isnormal( me )
            << " is0=" << ( me == 0 )
            << " is1=" << ( me == 1 )
            << " abs(ME)=" << std::abs( me )
            << " isnan=" << std::isnan( std::abs( me ) )
            << std::endl;
}

//============================================================================

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  // A base class encapsulating the calculation of event statistics on a CPU host or on a GPU device
  class CrossSectionKernelBase //: virtual public ICrossSectionKernel
  {
  protected:

    // Constructor from existing input and output buffers
    CrossSectionKernelBase( const BufferWeights& samplingWeights,       // input: sampling weights
                            const BufferMatrixElements& matrixElements, // input: matrix elements
                            EventStatistics& stats )                    // output: event statistics
      : m_samplingWeights( samplingWeights )
      , m_matrixElements( matrixElements )
      , m_stats( stats )
      , m_iter( 0 )
    {
      // NB: do not initialise EventStatistics (you may be asked to update an existing result)
    }

  public:

    // Destructor
    virtual ~CrossSectionKernelBase(){}

    // Update event statistics
    virtual void updateEventStatistics( const bool debug=false ) = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The buffer for the sampling weights
    const BufferWeights& m_samplingWeights;

    // The buffer for the output matrix elements
    const BufferMatrixElements& m_matrixElements;

    // The event statistics
    EventStatistics& m_stats;

    // The number of iterations processed so far
    size_t m_iter;

  };

  //--------------------------------------------------------------------------

  // A class encapsulating the calculation of event statistics on a CPU host
  class CrossSectionKernelHost final : public CrossSectionKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    CrossSectionKernelHost( const BufferWeights& samplingWeights,       // input: sampling weights
                            const BufferMatrixElements& matrixElements, // input: matrix elements
                            EventStatistics& stats,                     // output: event statistics
                            const size_t nevt );

    // Destructor
    virtual ~CrossSectionKernelHost(){}

    // Update event statistics
    void updateEventStatistics( const bool debug=false ) override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }

  };

  //--------------------------------------------------------------------------

  /*
#ifdef __CUDACC__
  // A class encapsulating the calculation of event statistics on a GPU device
  class CrossSectionKernelDevice : public CrossSectionKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    CrossSectionKernelDevice( const BufferWeights& samplingWeights,       // input: sampling weights
                              const BufferMatrixElements& matrixElements, // input: matrix elements
                              EventStatistics& stats,                     // output: event statistics
                              const size_t gpublocks,
                              const size_t gputhreads );

    // Destructor
    virtual ~CrossSectionKernelDevice(){}

    // Reset gpublocks and gputhreads
    void setGrid( const size_t gpublocks, const size_t gputhreads );

    // Update event statistics
    void updateEventStatistics( const bool debug=false ) override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;

  };
#endif
  */

  //--------------------------------------------------------------------------

}
#endif // CROSSSECTIONKERNELS_H
