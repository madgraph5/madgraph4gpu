#ifndef MemoryAccessCouplingsFixed_H
#define MemoryAccessCouplingsFixed_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

//#include "MemoryAccessHelpers.h"

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessCouplingsFixed
{
public:
  // TRIVIAL ACCESS to fixed-couplings buffers
  static __host__ __device__ inline cxtype*
  kernelAccess( fptype* buffer )
  {
    // FIXME: this assumes that ANY cxtype implementation is two adjacent fptypes - is it safer to cast to cxsmpl<fptype>*?
    return reinterpret_cast<cxtype*>( buffer );
  }

  // TRIVIAL ACCESS to fixed-couplings buffers
  static __host__ __device__ inline const cxtype*
  kernelAccessConst( const fptype* buffer )
  {
    // FIXME: this assumes that ANY cxtype implementation is two adjacent fptypes - is it safer to cast to cxsmpl<fptype>*?
    return reinterpret_cast<const cxtype*>( buffer );
  }
};

//----------------------------------------------------------------------------

typedef KernelAccessCouplingsFixed<false> HostAccessCouplingsFixed;
typedef KernelAccessCouplingsFixed<true> DeviceAccessCouplingsFixed;

//----------------------------------------------------------------------------

#endif // MemoryAccessCouplingsFixed_H
