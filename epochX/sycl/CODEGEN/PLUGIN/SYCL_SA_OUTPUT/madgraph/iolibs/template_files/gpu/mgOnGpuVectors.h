#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include "mgOnGpuTypes.h"
#include <iostream>

//==========================================================================

//------------------------------
// Vector types
//------------------------------
#if MGONGPU_MARRAY_DIM > 1
    #if defined MGONGPU_COMPLEX_CXSMPL
        typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
        typedef mgOnGpu::cxsmpl<fptype_sv> cxtype_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> int_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> bool_sv;
    #elif defined MGONGPU_COMPLEX_EXTRAS
        typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
        typedef extras::complex<fptype_sv> cxtype_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> int_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> bool_sv;
    #elif defined MGONGPU_COMPLEX_STD
        typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
        typedef std::complex<fptype_sv> cxtype_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> int_sv;
        typedef sycl::vec<long, MGONGPU_MARRAY_DIM> bool_sv;
    #else
        #error Unconfigured vector complex type. Add details to `mgOnGpuVectors.h` or set MGONGPU_MARRAY_DIM to 1 in `mgOnGpuConfig.h`.
    #endif
#else
    typedef fptype fptype_sv;
    typedef cxtype cxtype_sv;
    typedef long int_sv;
    typedef long bool_sv;
#endif

struct vector4 {
    fptype_sv w;
    fptype_sv x;
    fptype_sv y;
    fptype_sv z;
};

//------------------------------
// Vector operators
//------------------------------

#if MGONGPU_MARRAY_DIM > 1
    SYCL_EXTERNAL inline cxtype_sv operator*(const cxtype_sv& __x, const cxtype& __y) {
        cxtype_sv __r = __x;
        __r *= __y;
        return __r;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator*(const cxtype& __x, const cxtype_sv& __y) {
        return __y*__x;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator/(const cxtype_sv& __x, const cxtype& __y) {
        return __x/cxtype_sv(__y);
    }
    
    SYCL_EXTERNAL inline cxtype_sv& operator/=(cxtype_sv& __x, const cxtype& __y) {
        __x /= cxtype_sv(__y);
        return __x;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator/(const cxtype& __x, const cxtype_sv& __y) {
        return cxtype_sv(__x)/__y;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator+(const cxtype_sv& __x, const cxtype& __y) {
        cxtype_sv __r = __x;
        __r += __y;
        return __r;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator+(const cxtype& __x, const cxtype_sv& __y) {
        return __y + __x;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator-(const cxtype_sv& __x, const cxtype& __y) {
        cxtype_sv __r = __x;
        __r -= __y;
        return __r;
    }
    
    SYCL_EXTERNAL inline cxtype_sv operator-(const cxtype& __x, const cxtype_sv& __y) {
        return -__y + __x;
    }
#endif

//--------------------------------------------------------------------------


//==========================================================================
#endif // MGONGPUVECTORS_H
