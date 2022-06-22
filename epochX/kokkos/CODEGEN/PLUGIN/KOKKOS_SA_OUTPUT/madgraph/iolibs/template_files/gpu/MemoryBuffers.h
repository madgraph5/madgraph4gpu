#ifndef MemoryBuffers_H
#define MemoryBuffers_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuCxtypes.h"
#include "Kokkos_Core.hpp"

#include <sstream>

namespace mg5amcGpu
{
  //--------------------------------------------------------------------------

  // TEMPORARY? Take this from a PhysicsProcess class? Define them here directly in codegen?
  namespace MemoryBuffers
  {
    static constexpr size_t np4 = mgOnGpu::np4;
    static constexpr size_t nparf = mgOnGpu::nparf;
    static constexpr size_t npar = mgOnGpu::npar;
    static constexpr size_t nw6 = mgOnGpu::nw6;
    static constexpr size_t nx2 = mgOnGpu::nx2;
  }

  //--------------------------------------------------------------------------

  // An abstract interface encapsulating a given number of events
  class INumberOfEvents
  {
  public:
    virtual ~INumberOfEvents(){}
    virtual size_t nevt() const = 0;
  };

  //--------------------------------------------------------------------------

  // A class encapsulating a given number of events
  class NumberOfEvents : virtual public INumberOfEvents
  {
  public:
    NumberOfEvents( const size_t nevt )
      : m_nevt( nevt ){}
    virtual ~NumberOfEvents(){}
    virtual size_t nevt() const override { return m_nevt; }
  private:
    const size_t m_nevt;
  };

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer (not necessarily an event buffer)
  template<typename T>
  class BufferBase : virtual public INumberOfEvents
  {
  protected:
    BufferBase( const size_t size, const bool onDevice, sycl::queue q) : m_size( size ), m_data( nullptr ), m_isOnDevice( onDevice ), m_q( q ){}
    virtual ~BufferBase(){}
  public:
    T* data(){ return m_data; }
    const T* data() const{ return m_data; }
    sycl::queue get_queue(){ return m_q; }
    const sycl::queue get_queue() const{ return m_q; }
    T& operator[]( const size_t index ){ return m_data[index]; }
    const T& operator[]( const size_t index ) const { return m_data[index]; }
    size_t size() const{ return m_size; }
    size_t bytes() const{ return m_size * sizeof(T); }
    bool isOnDevice() const { return m_isOnDevice; }
    virtual size_t nevt() const override { throw std::runtime_error( "This BufferBase is not an event buffer" ); }
  protected:
    const size_t m_size;
    T* m_data;
    const bool m_isOnDevice;
    sycl::queue m_q;
 };

  //--------------------------------------------------------------------------

  // A class encapsulating a SYCL pinned host buffer
  template<typename T>
  class PinnedHostBufferBase : public BufferBase<T>
  {
  public:
    PinnedHostBufferBase( const size_t size, sycl::queue q ) : BufferBase<T>( size, false, q )
    {
      this->m_data = sycl::malloc_host<T>(this->size(), this->m_q);
    }
    virtual ~PinnedHostBufferBase()
    {
      sycl::free( this->m_data, this->m_q );
    }
  };

  //--------------------------------------------------------------------------

  // A class encapsulating a SYCL device buffer
  template<typename T>
  class DeviceBufferBase : public BufferBase<T>
  {
  public:
    DeviceBufferBase( const size_t size, sycl::queue q ) : BufferBase<T>( size, true, q )
    {
      this->m_data = sycl::malloc_device<T>(this->size(), this->m_q);
    }
    virtual ~DeviceBufferBase()
    {
      sycl::free( this->m_data, this->m_q );
    }
  };

  //--------------------------------------------------------------------------

  // A class encapsulating a SYCL pinned host buffer for a given number of events
  template<typename T, size_t sizePerEvent>
  class PinnedHostBuffer : public PinnedHostBufferBase<T>, virtual private NumberOfEvents
  {
  public:
    PinnedHostBuffer( const size_t nevt, sycl::queue q)
      : NumberOfEvents( nevt )
      , PinnedHostBufferBase<T>( sizePerEvent * nevt, q ){}
    virtual ~PinnedHostBuffer(){}
    virtual size_t nevt() const override final { return NumberOfEvents::nevt(); }
  };

  //--------------------------------------------------------------------------

  // A class encapsulating a SYCL device buffer for a given number of events
  template<typename T, size_t sizePerEvent>
  class DeviceBuffer : public DeviceBufferBase<T>, virtual private NumberOfEvents
  {
  public:
    DeviceBuffer( const size_t nevt, sycl::queue q )
      : NumberOfEvents( nevt )
      , DeviceBufferBase<T>( sizePerEvent * nevt, q ){}
    virtual ~DeviceBuffer(){}
    virtual size_t nevt() const override final { return NumberOfEvents::nevt(); }
  };

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for random numbers
  typedef BufferBase<fptype> BufferRandomNumbers;

  // The size (number of elements) per event in a memory buffer for random numbers
  constexpr size_t sizePerEventRandomNumbers = MemoryBuffers::np4 * MemoryBuffers::nparf;

  // A class encapsulating a SYCL pinned host buffer for random numbers
  typedef PinnedHostBuffer<fptype, sizePerEventRandomNumbers> PinnedHostBufferRandomNumbers;
  // A class encapsulating a SYCL device buffer for random numbers
  typedef DeviceBuffer<fptype, sizePerEventRandomNumbers> DeviceBufferRandomNumbers;

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for momenta
  typedef BufferBase<fptype> BufferMomenta;

  // The size (number of elements) per event in a memory buffer for momenta
  constexpr size_t sizePerEventMomenta = MemoryBuffers::np4 * MemoryBuffers::npar;

  // A class encapsulating a SYCL pinned host buffer for momenta
  typedef PinnedHostBuffer<fptype, sizePerEventMomenta> PinnedHostBufferMomenta;
  // A class encapsulating a SYCL device buffer for momenta
  typedef DeviceBuffer<fptype, sizePerEventMomenta> DeviceBufferMomenta;

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for sampling weights
  typedef BufferBase<fptype> BufferWeights;

  // The size (number of elements) per event in a memory buffer for sampling weights
  constexpr size_t sizePerEventWeights = 1;

  // A class encapsulating a SYCL pinned host buffer for sampling weights
  typedef PinnedHostBuffer<fptype, sizePerEventWeights> PinnedHostBufferWeights;
  // A class encapsulating a SYCL device buffer for sampling weights
  typedef DeviceBuffer<fptype, sizePerEventWeights> DeviceBufferWeights;

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for matrix elements
  typedef BufferBase<fptype> BufferMatrixElements;

  // The size (number of elements) per event in a memory buffer for matrix elements
  constexpr size_t sizePerEventMatrixElements = 1;

  // A class encapsulating a SYCL pinned host buffer for matrix elements
  typedef PinnedHostBuffer<fptype, sizePerEventMatrixElements> PinnedHostBufferMatrixElements;
  // A class encapsulating a SYCL device buffer for matrix elements
  typedef DeviceBuffer<fptype, sizePerEventMatrixElements> DeviceBufferMatrixElements;

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for the helicity mask
  typedef BufferBase<bool> BufferHelicityMask;

  // A class encapsulating a SYCL pinned host buffer for the helicity mask
  typedef PinnedHostBufferBase<bool> PinnedHostBufferHelicityMask;
  // A class encapsulating a SYCL device buffer for the helicity mask
  typedef DeviceBufferBase<bool> DeviceBufferHelicityMask;

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for wavefunctions
  typedef BufferBase<fptype> BufferWavefunctions;

  // The size (number of elements) per event in a memory buffer for wavefunctions
  constexpr size_t sizePerEventWavefunctions = MemoryBuffers::nw6 * MemoryBuffers::nx2;

  // A class encapsulating a SYCL pinned host buffer for wavefunctions
  typedef PinnedHostBuffer<fptype, sizePerEventWavefunctions> PinnedHostBufferWavefunctions;
  // A class encapsulating a SYCL device buffer for wavefunctions
  typedef DeviceBuffer<fptype, sizePerEventWavefunctions> DeviceBufferWavefunctions;

  //--------------------------------------------------------------------------

  template<class Tdst, class Tsrc>
  void copyDeviceFromHost( Tdst& dst, const Tsrc& src ) // keep the same order of arguments as in memcpy
  {
    if ( dst.size() != src.size() )
    {
      std::ostringstream sstr;
      sstr << "Size (#elements) mismatch in copyDeviceFromHost: dst=" << dst.size() << ", src=" << src.size();
      throw std::runtime_error( sstr.str() );
    }
    if ( dst.bytes() != src.bytes() )
    {
      std::ostringstream sstr;
      sstr << "Size (#bytes) mismatch in copyDeviceFromHost: dst=" << dst.bytes() << ", src=" << src.bytes();
      throw std::runtime_error( sstr.str() );
    }
    sycl::queue q = dst.get_queue();
    q.memcpy( dst.data(), src.data(), src.bytes() ).wait();
  }

  //--------------------------------------------------------------------------

  template<class Tdst, class Tsrc>
  void copyHostFromDevice( Tdst& dst, const Tsrc& src ) // keep the same order of arguments as in memcpy
  {
    if ( dst.size() != src.size() )
    {
      std::ostringstream sstr;
      sstr << "Size (#elements) mismatch in copyHostFromDevice: dst=" << dst.size() << ", src=" << src.size();
      throw std::runtime_error( sstr.str() );
    }
    if ( dst.bytes() != src.bytes() )
    {
      std::ostringstream sstr;
      sstr << "Size (#bytes) mismatch in copyHostFromDevice: dst=" << dst.bytes() << ", src=" << src.bytes();
      throw std::runtime_error( sstr.str() );
    }
    sycl::queue q = src.get_queue();
    q.memcpy( dst.data(), src.data(), src.bytes() ).wait();
  }

  //--------------------------------------------------------------------------
}

#endif // MemoryBuffers_H
