#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#ifndef __CUDACC__

#include "mgOnGpuTypes.h"

namespace mgOnGpu
{
  const int neppV = neppM;

  // --- Type definitions
  typedef fptype fptype_v[neppV]; // RRRR
  typedef cxtype cxtype_v[neppV]; // RIRIRIRI (not RRRRIIII)

  // --- Type definitions (using vector compiler extensions)
  //typedef fptype fptype_v __attribute__ ((vector_size (neppV * sizeof(fptype)))); // OK
  //typedef cxtype cxtype_v __attribute__ ((vector_size (neppV * sizeof(cxtype)))); // invalid vector type
}

// Expose typedefs and operators outside the namespace
using mgOnGpu::neppV;
using mgOnGpu::fptype_v;
using mgOnGpu::cxtype_v;

#endif

#endif // MGONGPUVECTORS_H
