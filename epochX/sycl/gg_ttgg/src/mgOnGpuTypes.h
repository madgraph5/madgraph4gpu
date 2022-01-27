#ifndef MGONGPUTYPES_H
#define MGONGPUTYPES_H 1

#include "mgOnGpuConfig.h"
#include "extras.h"
#include <cmath>
#include <complex>

namespace mgOnGpu
{

// --- Type definitions

// Complex type: cxtype
typedef extras::complex<fptype> cxtype; // two doubles: RI

}

// Expose typedefs and operators outside the namespace
using mgOnGpu::cxtype;

// --- Functions and operators for floating point types

inline const fptype& fpmax( const fptype& a, const fptype& b ) {
    return std::max( a, b );
}

inline const fptype& fpmin( const fptype& a, const fptype& b ) {
    return std::min( a, b );
}

inline fptype fpsqrt( const fptype& f ) {
    return sycl::sqrt( f );
}

// --- Functions and operators for complex types

inline cxtype cxmake( const fptype& r, const fptype& i ) {
    return cxtype( r, i ); // std::complex<fptype> constructor
}

inline fptype cxreal( const cxtype& c ) {
    return c.real(); // std::complex<fptype>::real()
}

inline fptype cximag( const cxtype& c ) {
    return c.imag(); // std::complex<fptype>::imag()
}

inline const cxtype& cxmake( const cxtype& c ) { // std::complex to std::complex (float-to-float or double-to-double)
    return c;
}

inline cxtype cxmake( const std::complex<double>& c ) { // std::complex to std::complex (cast double-to-fptype)
    return cxmake( (fptype)c.real(), (fptype)c.imag() );
}


// SYCL complex operator overloads
inline std::ostream& operator<<( std::ostream& out, const cxtype& c ) {
    out << "[" << cxreal(c) << "," << cximag(c) << "]";
    return out;
}

inline SYCL_EXTERNAL cxtype operator+( const cxtype a ) {
    return a;
}

inline SYCL_EXTERNAL cxtype operator-( const cxtype& a ) {
    return cxmake( -cxreal(a), -cximag(a) );
}

inline SYCL_EXTERNAL cxtype operator+( const fptype& a, const cxtype& b ) {
    return cxmake( a, 0 ) + b;
}

inline SYCL_EXTERNAL cxtype operator-( const fptype& a, const cxtype& b ) {
    return cxmake( a, 0 ) - b;
}

inline SYCL_EXTERNAL cxtype operator*( const fptype& a, const cxtype& b ) {
    return cxmake( a, 0 ) * b;
}

inline SYCL_EXTERNAL cxtype conj( const cxtype& c ) {
    return cxmake( cxreal( c ), -cximag( c ) );
}

inline SYCL_EXTERNAL cxtype cxconj( const cxtype& c ) {
    return conj( c ); 
}
#endif // MGONGPUTYPES_H
