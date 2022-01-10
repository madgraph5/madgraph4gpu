#ifndef MemoryBuffers_H
#define MemoryBuffers_H 1

#include "mgOnGpuConfig.h"

#include "checkCuda.h"

#include <sstream>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer
  template<typename T>
  class BufferBase
  {
  protected:
    BufferBase( const size_t size ) : m_size( size ), m_data( nullptr ){}
    virtual ~BufferBase(){}
  public:
    T* data(){ return m_data; }
    const T* data() const{ return m_data; }
    size_t size() const{ return m_size; }
    size_t bytes() const{ return m_size * sizeof(T); }
  protected:
    const size_t m_size;
    T* m_data;
  };

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A base class encapsulating a C++ host buffer
  template<typename T>
  class HostBufferBase : public BufferBase<T>
  {
  protected:
    HostBufferBase( const size_t size ) : BufferBase<T>( size )
    {
      this->m_data = new( std::align_val_t{ cppAlign } ) T[ size ]();
    }    
    virtual ~HostBufferBase()
    {
      ::operator delete( this->m_data, std::align_val_t{ cppAlign } );
    }    
  public:
    static constexpr int cppAlign = mgOnGpu::cppAlign;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A base class encapsulating a CUDA pinned host buffer
  template<typename T>
  class PinnedHostBufferBase : public BufferBase<T>
  {
  protected:
    PinnedHostBufferBase( const size_t size ) : BufferBase<T>( size )
    {
      checkCuda( cudaMallocHost( &(this->m_data), this->bytes() ) );
    }    
    virtual ~PinnedHostBufferBase()
    {
      checkCuda( cudaFreeHost( this->m_data ) );
    }    
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A base class encapsulating a CUDA device buffer
  template<typename T>
  class DeviceBufferBase : public BufferBase<T>
  {
  protected:
    DeviceBufferBase( const size_t size ) : BufferBase<T>( size )
    {
      checkCuda( cudaMalloc( &(this->m_data), this->bytes() ) );
    }    
    virtual ~DeviceBufferBase()
    {
      checkCuda( cudaFree( this->m_data ) );
    }    
  };
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for random numbers
  typedef BufferBase<fptype> BufferRandomNumbers;

  //--------------------------------------------------------------------------

  // TEMPORARY? Take this from a PhysicsProcess class? Define them here directly in codegen?
  namespace MemoryBuffers
  {
    static constexpr int nparf = mgOnGpu::nparf;
    static constexpr int np4 = mgOnGpu::np4;
  }
  
  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for random numbers
  class HostBufferRandomNumbers : public HostBufferBase<fptype>
  {
  public:
    HostBufferRandomNumbers( const int nevt )
      : HostBufferBase<fptype>( MemoryBuffers::np4 * MemoryBuffers::nparf * nevt )
      , m_nevt( nevt ){}    
    virtual ~HostBufferRandomNumbers(){}
    size_t nevt() const{ return m_nevt; }
  private:
    const size_t m_nevt;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA pinned host buffer for random numbers
  class PinnedHostBufferRandomNumbers : public PinnedHostBufferBase<fptype>
  {
  public:
    PinnedHostBufferRandomNumbers( const int nevt )
      : PinnedHostBufferBase<fptype>( MemoryBuffers::np4 * MemoryBuffers::nparf * nevt )
      , m_nevt( nevt ){}    
    virtual ~PinnedHostBufferRandomNumbers(){}
    size_t nevt() const{ return m_nevt; }
  private:
    const size_t m_nevt;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA device buffer for random numbers
  class DeviceBufferRandomNumbers : public DeviceBufferBase<fptype>
  {
  public:
    DeviceBufferRandomNumbers( const int nevt )
      : DeviceBufferBase<fptype>( MemoryBuffers::np4 * MemoryBuffers::nparf * nevt )
      , m_nevt( nevt ){}    
    virtual ~DeviceBufferRandomNumbers(){}
    size_t nevt() const{ return m_nevt; }
  private:
    const size_t m_nevt;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
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
    // NB (PR #45): cudaMemcpy involves an intermediate memcpy to pinned memory if host array is a not a pinned host array
    checkCuda( cudaMemcpy( dst.data(), src.data(), src.bytes(), cudaMemcpyHostToDevice ) );
  }
#endif

  //--------------------------------------------------------------------------
  
#ifdef __CUDACC__
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
    // NB (PR #45): cudaMemcpy involves an intermediate memcpy to pinned memory if host array is a not a pinned host array
    checkCuda( cudaMemcpy( dst.data(), src.data(), src.bytes(), cudaMemcpyDeviceToHost ) );
  }
#endif

  //--------------------------------------------------------------------------
}

#endif // MemoryBuffers_H
