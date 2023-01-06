#ifndef MGONGPUVECTORS_H
#define MGONGPUVECTORS_H 1

#include "mgOnGpuTypes.h"
#include <iostream>

//==========================================================================

//------------------------------
// Vector types
//------------------------------
#if defined MGONGPU_COMPLEX_CXSMPLVEC
    typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
    typedef mgOnGpu::cxsmpl<fptype_sv> cxtype_sv;
#elif defined MGONGPU_COMPLEX_EXTRASVEC
    typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
    typedef extras::complex<fptype_sv> cxtype_sv;
#elif defined MGONGPU_COMPLEX_STDVEC
    typedef sycl::vec<fptype, MGONGPU_MARRAY_DIM> fptype_sv;
    typedef std::complex<fptype_sv> cxtype_sv;
#else
    typedef fptype fptype_sv;
    typedef cxtype cxtype_sv;
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

#if defined MGONGPU_USE_VEC
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
