#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include <iostream>

//------------------------------
// Vector types
//------------------------------

#include "mgOnGpuTypes.h"

namespace mgOnGpu {
#ifdef MGONGPU_NEPPV
    const int neppV = MGONGPU_NEPPV;
#else
    const int neppV = 1;
#endif
}

// Expose typedefs outside the namespace
using mgOnGpu::neppV;

// Printout to std::cout for user defined types
inline void print( const fptype& f ) { std::cout << f << std::endl; }

inline void print( const cxtype& c ) { std::cout << c << std::endl; }

// Operators for bool
inline fptype fpternary( const bool& mask, const fptype& a, const fptype& b ) {
    return ( mask ? a : b );
}

inline cxtype cxternary( const bool& mask, const cxtype& a, const cxtype& b ) {
    return ( mask ? a : b );
}

inline bool maskor( const bool& mask ) {
  return mask;
}

//--------------------------------------------------------------------------

// Scalar-or-vector types:
typedef bool bool_sv;
typedef fptype fptype_sv;
typedef cxtype cxtype_sv;

// Scalar-or-vector zeros:
inline cxtype cxzero_sv(){ return cxtype{ fptype{0}, fptype{0} }; }

//--------------------------------------------------------------------------

#endif // MGONGPUVECTORS_H
