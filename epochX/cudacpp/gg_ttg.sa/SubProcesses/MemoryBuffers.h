// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2021, based on earlier work by S. Hageboeck) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Roiser, A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

#ifndef MemoryBuffers_H
#define MemoryBuffers_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "CudaRuntime.h"
#include "Parameters_sm.h"

#include <sstream>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
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
    static constexpr size_t ndcoup = Parameters_sm_dependentCouplings::ndcoup;
  }

  //--------------------------------------------------------------------------

  // An abstract interface encapsulating a given number of events
  class INumberOfEvents
  {
  public:
    virtual ~INumberOfEvents() {}
    virtual size_t nevt() const = 0;
  };

  //--------------------------------------------------------------------------

  // A class encapsulating a given number of events
  class NumberOfEvents : virtual public INumberOfEvents
  {
  public:
    NumberOfEvents( const size_t nevt )
      : m_nevt( nevt ) {}
    virtual ~NumberOfEvents() {}
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
    BufferBase( const size_t size, const bool onDevice )
      : m_size( size ), m_data( nullptr ), m_isOnDevice( onDevice ) {}
    virtual ~BufferBase() {}
  public:
    T* data() { return m_data; }
    const T* data() const { return m_data; }
    T& operator[]( const size_t index ) { return m_data[index]; }
    const T& operator[]( const size_t index ) const { return m_data[index]; }
    size_t size() const { return m_size; }
    size_t bytes() const { return m_size * sizeof( T ); }
    bool isOnDevice() const { return m_isOnDevice; }
    virtual size_t nevt() const override { throw std::runtime_error( "This BufferBase is not an event buffer" ); }
  protected:
    const size_t m_size;
    T* m_data;
    const bool m_isOnDevice;
  };

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  constexpr bool HostBufferALIGNED = false;   // ismisaligned=false
  constexpr bool HostBufferMISALIGNED = true; // ismisaligned=true

  // A class encapsulating a C++ host buffer
  template<typename T, bool ismisaligned>
  class HostBufferBase : public BufferBase<T>
  {
  public:
    HostBufferBase( const size_t size )
      : BufferBase<T>( size, false )
    {
      if constexpr( !ismisaligned )
        this->m_data = new( std::align_val_t( cppAlign ) ) T[size]();
      else
        this->m_data = new( std::align_val_t( cppAlign ) ) T[size + 1]() + 1; // TEST MISALIGNMENT!
    }
    virtual ~HostBufferBase()
    {
      if constexpr( !ismisaligned )
        ::operator delete[]( this->m_data, std::align_val_t( cppAlign ) );
      else
        ::operator delete[]( ( this->m_data ) - 1, std::align_val_t( cppAlign ) ); // TEST MISALIGNMENT!
    }
    static constexpr bool isaligned() { return !ismisaligned; }
  public:
    static constexpr size_t cppAlign = mgOnGpu::cppAlign;
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA pinned host buffer
  template<typename T>
  class PinnedHostBufferBase : public BufferBase<T>
  {
  public:
    PinnedHostBufferBase( const size_t size )
      : BufferBase<T>( size, false )
    {
      checkCuda( cudaMallocHost( &( this->m_data ), this->bytes() ) );
    }
    virtual ~PinnedHostBufferBase()
    {
      checkCuda( cudaFreeHost( this->m_data ) );
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA device buffer
  template<typename T>
  class DeviceBufferBase : public BufferBase<T>
  {
  public:
    DeviceBufferBase( const size_t size )
      : BufferBase<T>( size, true )
    {
      checkCuda( cudaMalloc( &( this->m_data ), this->bytes() ) );
    }
    virtual ~DeviceBufferBase()
    {
      checkCuda( cudaFree( this->m_data ) );
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for a given number of events
  template<typename T, size_t sizePerEvent, bool ismisaligned>
  class HostBuffer : public HostBufferBase<T, ismisaligned>, virtual private NumberOfEvents
  {
  public:
    HostBuffer( const size_t nevt )
      : NumberOfEvents( nevt )
      , HostBufferBase<T, ismisaligned>( sizePerEvent * nevt ) {}
    virtual ~HostBuffer() {}
    virtual size_t nevt() const override final { return NumberOfEvents::nevt(); }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA pinned host buffer for a given number of events
  template<typename T, size_t sizePerEvent>
  class PinnedHostBuffer : public PinnedHostBufferBase<T>, virtual private NumberOfEvents
  {
  public:
    PinnedHostBuffer( const size_t nevt )
      : NumberOfEvents( nevt )
      , PinnedHostBufferBase<T>( sizePerEvent * nevt ) {}
    virtual ~PinnedHostBuffer() {}
    virtual size_t nevt() const override final { return NumberOfEvents::nevt(); }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  // A class encapsulating a CUDA device buffer for a given number of events
  template<typename T, size_t sizePerEvent>
  class DeviceBuffer : public DeviceBufferBase<T>, virtual private NumberOfEvents
  {
  public:
    DeviceBuffer( const size_t nevt )
      : NumberOfEvents( nevt )
      , DeviceBufferBase<T>( sizePerEvent * nevt ) {}
    virtual ~DeviceBuffer() {}
    virtual size_t nevt() const override final { return NumberOfEvents::nevt(); }
  };
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for momenta random numbers
  typedef BufferBase<fptype> BufferRndNumMomenta;

  // The size (number of elements) per event in a memory buffer for momenta random numbers
  constexpr size_t sizePerEventRndNumMomenta = MemoryBuffers::np4 * MemoryBuffers::nparf;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for momenta random numbers
  typedef HostBuffer<fptype, sizePerEventRndNumMomenta, HostBufferALIGNED> HostBufferRndNumMomenta;
#else
  // A class encapsulating a CUDA pinned host buffer for momenta random numbers
  typedef PinnedHostBuffer<fptype, sizePerEventRndNumMomenta> PinnedHostBufferRndNumMomenta;
  // A class encapsulating a CUDA device buffer for momenta random numbers
  typedef DeviceBuffer<fptype, sizePerEventRndNumMomenta> DeviceBufferRndNumMomenta;
#endif

  //--------------------------------------------------------------------------

  /*
  // A base class encapsulating a memory buffer with ONE fptype per event
  typedef BufferBase<fptype> BufferOneFp;

  // The size (number of elements) per event in a memory buffer with ONE fptype per event
  constexpr size_t sizePerEventOneFp = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer with ONE fptype per event
  typedef HostBuffer<fptype, sizePerEventOneFp, HostBufferALIGNED> HostBufferOneFp;
#else
  // A class encapsulating a CUDA pinned host buffer for gs
  typedef PinnedHostBuffer<fptype, sizePerEventOneFp> PinnedHostBufferOneFp;
  // A class encapsulating a CUDA device buffer for gs
  typedef DeviceBuffer<fptype, sizePerEventOneFp> DeviceBufferOneFp;
#endif

  // Memory buffers for Gs (related to the event-by-event strength of running coupling constant alphas QCD)
  typedef BufferOneFp BufferGs;
  typedef HostBufferOneFp HostBufferGs;
  typedef PinnedHostBufferOneFp PinnedHostBufferGs;
  typedef DeviceBufferOneFp DeviceBufferGs;
  */

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for Gs (related to the event-by-event strength of running coupling constant alphas QCD)
  typedef BufferBase<fptype> BufferGs;

  // The size (number of elements) per event in a memory buffer for Gs
  constexpr size_t sizePerEventGs = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for gs
  typedef HostBuffer<fptype, sizePerEventGs, HostBufferALIGNED> HostBufferGs;
#else
  // A class encapsulating a CUDA pinned host buffer for gs
  typedef PinnedHostBuffer<fptype, sizePerEventGs> PinnedHostBufferGs;
  // A class encapsulating a CUDA device buffer for gs
  typedef DeviceBuffer<fptype, sizePerEventGs> DeviceBufferGs;
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  // A base class encapsulating a memory buffer for numerators (of the multichannel single-diagram enhancement factors)
  typedef BufferBase<fptype> BufferNumerators;

  // The size (number of elements) per event in a memory buffer for numerators
  constexpr size_t sizePerEventNumerators = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for gs
  typedef HostBuffer<fptype, sizePerEventNumerators, HostBufferALIGNED> HostBufferNumerators;
#else
  // A class encapsulating a CUDA pinned host buffer for gs
  typedef PinnedHostBuffer<fptype, sizePerEventNumerators> PinnedHostBufferNumerators;
  // A class encapsulating a CUDA device buffer for gs
  typedef DeviceBuffer<fptype, sizePerEventNumerators> DeviceBufferNumerators;
#endif
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  // A base class encapsulating a memory buffer for denominators (of the multichannel single-diagram enhancement factors)
  typedef BufferBase<fptype> BufferDenominators;

  // The size (number of elements) per event in a memory buffer for denominators
  constexpr size_t sizePerEventDenominators = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for gs
  typedef HostBuffer<fptype, sizePerEventDenominators, HostBufferALIGNED> HostBufferDenominators;
#else
  // A class encapsulating a CUDA pinned host buffer for gs
  typedef PinnedHostBuffer<fptype, sizePerEventDenominators> PinnedHostBufferDenominators;
  // A class encapsulating a CUDA device buffer for gs
  typedef DeviceBuffer<fptype, sizePerEventDenominators> DeviceBufferDenominators;
#endif
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for couplings that depend on the event-by-event running coupling constant alphas QCD
  typedef BufferBase<fptype> BufferCouplings;

  // The size (number of elements) per event in a memory buffer for random numbers
  constexpr size_t sizePerEventCouplings = MemoryBuffers::ndcoup * MemoryBuffers::nx2;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for gs
  typedef HostBuffer<fptype, sizePerEventCouplings, HostBufferALIGNED> HostBufferCouplings;
#else
  // A class encapsulating a CUDA pinned host buffer for gs
  typedef PinnedHostBuffer<fptype, sizePerEventCouplings> PinnedHostBufferCouplings;
  // A class encapsulating a CUDA device buffer for gs
  typedef DeviceBuffer<fptype, sizePerEventCouplings> DeviceBufferCouplings;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for momenta
  typedef BufferBase<fptype> BufferMomenta;

  // The size (number of elements) per event in a memory buffer for momenta
  constexpr size_t sizePerEventMomenta = MemoryBuffers::np4 * MemoryBuffers::npar;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for momenta
  typedef HostBuffer<fptype, sizePerEventMomenta, HostBufferALIGNED> HostBufferMomenta;
  //typedef HostBuffer<fptype, sizePerEventMomenta, HostBufferMISALIGNED> HostBufferMomenta; // TEST MISALIGNMENT!
#else
  // A class encapsulating a CUDA pinned host buffer for momenta
  typedef PinnedHostBuffer<fptype, sizePerEventMomenta> PinnedHostBufferMomenta;
  // A class encapsulating a CUDA device buffer for momenta
  typedef DeviceBuffer<fptype, sizePerEventMomenta> DeviceBufferMomenta;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for sampling weights
  typedef BufferBase<fptype> BufferWeights;

  // The size (number of elements) per event in a memory buffer for sampling weights
  constexpr size_t sizePerEventWeights = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for sampling weights
  typedef HostBuffer<fptype, sizePerEventWeights, HostBufferALIGNED> HostBufferWeights;
#else
  // A class encapsulating a CUDA pinned host buffer for sampling weights
  typedef PinnedHostBuffer<fptype, sizePerEventWeights> PinnedHostBufferWeights;
  // A class encapsulating a CUDA device buffer for sampling weights
  typedef DeviceBuffer<fptype, sizePerEventWeights> DeviceBufferWeights;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for matrix elements
  typedef BufferBase<fptype> BufferMatrixElements;

  // The size (number of elements) per event in a memory buffer for matrix elements
  constexpr size_t sizePerEventMatrixElements = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for matrix elements
  typedef HostBuffer<fptype, sizePerEventMatrixElements, HostBufferALIGNED> HostBufferMatrixElements;
#else
  // A class encapsulating a CUDA pinned host buffer for matrix elements
  typedef PinnedHostBuffer<fptype, sizePerEventMatrixElements> PinnedHostBufferMatrixElements;
  // A class encapsulating a CUDA device buffer for matrix elements
  typedef DeviceBuffer<fptype, sizePerEventMatrixElements> DeviceBufferMatrixElements;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for the helicity mask
  typedef BufferBase<bool> BufferHelicityMask;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for the helicity mask
  typedef HostBufferBase<bool, HostBufferALIGNED> HostBufferHelicityMask;
#else
  // A class encapsulating a CUDA pinned host buffer for the helicity mask
  typedef PinnedHostBufferBase<bool> PinnedHostBufferHelicityMask;
  // A class encapsulating a CUDA device buffer for the helicity mask
  typedef DeviceBufferBase<bool> DeviceBufferHelicityMask;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for wavefunctions
  typedef BufferBase<fptype> BufferWavefunctions;

  // The size (number of elements) per event in a memory buffer for wavefunctions
  constexpr size_t sizePerEventWavefunctions = MemoryBuffers::nw6 * MemoryBuffers::nx2;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for wavefunctions
  typedef HostBuffer<fptype, sizePerEventWavefunctions, HostBufferALIGNED> HostBufferWavefunctions;
#else
  // A class encapsulating a CUDA pinned host buffer for wavefunctions
  typedef PinnedHostBuffer<fptype, sizePerEventWavefunctions> PinnedHostBufferWavefunctions;
  // A class encapsulating a CUDA device buffer for wavefunctions
  typedef DeviceBuffer<fptype, sizePerEventWavefunctions> DeviceBufferWavefunctions;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for helicity random numbers
  typedef BufferBase<fptype> BufferRndNumHelicity;

  // The size (number of elements) per event in a memory buffer for helicity random numbers
  constexpr size_t sizePerEventRndNumHelicity = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for helicity random numbers
  typedef HostBuffer<fptype, sizePerEventRndNumHelicity, HostBufferALIGNED> HostBufferRndNumHelicity;
#else
  // A class encapsulating a CUDA pinned host buffer for helicity random numbers
  typedef PinnedHostBuffer<fptype, sizePerEventRndNumHelicity> PinnedHostBufferRndNumHelicity;
  // A class encapsulating a CUDA device buffer for helicity random numbers
  typedef DeviceBuffer<fptype, sizePerEventRndNumHelicity> DeviceBufferRndNumHelicity;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for color random numbers
  typedef BufferBase<fptype> BufferRndNumColor;

  // The size (number of elements) per event in a memory buffer for color random numbers
  constexpr size_t sizePerEventRndNumColor = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for color random numbers
  typedef HostBuffer<fptype, sizePerEventRndNumColor, HostBufferALIGNED> HostBufferRndNumColor;
#else
  // A class encapsulating a CUDA pinned host buffer for color random numbers
  typedef PinnedHostBuffer<fptype, sizePerEventRndNumColor> PinnedHostBufferRndNumColor;
  // A class encapsulating a CUDA device buffer for color random numbers
  typedef DeviceBuffer<fptype, sizePerEventRndNumColor> DeviceBufferRndNumColor;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for helicity selection
  typedef BufferBase<int> BufferSelectedHelicity;

  // The size (number of elements) per event in a memory buffer for helicity selection
  constexpr size_t sizePerEventSelectedHelicity = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for helicity selection
  typedef HostBuffer<int, sizePerEventSelectedHelicity, HostBufferALIGNED> HostBufferSelectedHelicity;
#else
  // A class encapsulating a CUDA pinned host buffer for helicity selection
  typedef PinnedHostBuffer<int, sizePerEventSelectedHelicity> PinnedHostBufferSelectedHelicity;
  // A class encapsulating a CUDA device buffer for helicity selection
  typedef DeviceBuffer<int, sizePerEventSelectedHelicity> DeviceBufferSelectedHelicity;
#endif

  //--------------------------------------------------------------------------

  // A base class encapsulating a memory buffer for color selection
  typedef BufferBase<int> BufferSelectedColor;

  // The size (number of elements) per event in a memory buffer for color selection
  constexpr size_t sizePerEventSelectedColor = 1;

#ifndef __CUDACC__
  // A class encapsulating a C++ host buffer for color selection
  typedef HostBuffer<int, sizePerEventSelectedColor, HostBufferALIGNED> HostBufferSelectedColor;
#else
  // A class encapsulating a CUDA pinned host buffer for color selection
  typedef PinnedHostBuffer<int, sizePerEventSelectedColor> PinnedHostBufferSelectedColor;
  // A class encapsulating a CUDA device buffer for color selection
  typedef DeviceBuffer<int, sizePerEventSelectedColor> DeviceBufferSelectedColor;
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  template<class Tdst, class Tsrc>
  void copyDeviceFromHost( Tdst& dst, const Tsrc& src ) // keep the same order of arguments as in memcpy
  {
    if( dst.size() != src.size() )
    {
      std::ostringstream sstr;
      sstr << "Size (#elements) mismatch in copyDeviceFromHost: dst=" << dst.size() << ", src=" << src.size();
      throw std::runtime_error( sstr.str() );
    }
    if( dst.bytes() != src.bytes() )
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
    if( dst.size() != src.size() )
    {
      std::ostringstream sstr;
      sstr << "Size (#elements) mismatch in copyHostFromDevice: dst=" << dst.size() << ", src=" << src.size();
      throw std::runtime_error( sstr.str() );
    }
    if( dst.bytes() != src.bytes() )
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
