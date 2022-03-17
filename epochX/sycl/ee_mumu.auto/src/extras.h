#ifndef EXTRAS_H
#define EXTRAS_H 1

#include <CL/sycl.hpp>

namespace extras {

// Forward declarations
template<typename FPType> class complex;

/// Return complex conjugate
template<typename FPType>
SYCL_EXTERNAL complex<FPType> conj(const complex<FPType>&);

template< typename FPType >
class complex {
public:
  
    //Constructors
    SYCL_EXTERNAL constexpr complex(const FPType& __r = FPType(), const FPType& __i = FPType())
        : r_(__r), i_(__i) { }
    SYCL_EXTERNAL constexpr complex(const complex&) = default;
  
    //Converter
    template<typename CXType>
    SYCL_EXTERNAL constexpr complex(const complex<CXType>& __z)
        : r_(__z.real()), i_(__z.imag()) { }
  
    //Data Accessors
    SYCL_EXTERNAL constexpr FPType real() const { return r_; }
    SYCL_EXTERNAL constexpr FPType imag() const { return i_; }
  
  
    //Operators
    /// Assign scalar
    SYCL_EXTERNAL complex<FPType>& operator=(const FPType&);
  
    /// Add scalar
    SYCL_EXTERNAL complex<FPType>& operator+=(const FPType& __t) {
      r_ += __t;
      return *this;
    }
  
    /// Subtract scalar
    SYCL_EXTERNAL complex<FPType>& operator-=(const FPType& __t) {
      r_ -= __t;
      return *this;
    }
  
    /// Multiply by scalar
    SYCL_EXTERNAL complex<FPType>& operator*=(const FPType&);
  
    /// Divide by scalar
    SYCL_EXTERNAL complex<FPType>& operator/=(const FPType&);
  
    /// Compiler generated assignment operator
    SYCL_EXTERNAL complex& operator=(const complex&) = default;
  
    /// Assign another complex number to this one.
    template<typename CXType>
    SYCL_EXTERNAL complex<FPType>& operator=(const complex<CXType>&);
    
    /// Add complex number
    template<typename CXType>
    SYCL_EXTERNAL complex<FPType>& operator+=(const complex<CXType>&);
  
    /// Subtract complex number
    template<typename CXType>
    SYCL_EXTERNAL complex<FPType>& operator-=(const complex<CXType>&);
  
    /// Multiply by complex number
    template<typename CXType>
    SYCL_EXTERNAL complex<FPType>& operator*=(const complex<CXType>&);
  
    /// Divide by complex number
    template<typename CXType>
    SYCL_EXTERNAL complex<FPType>& operator/=(const complex<CXType>&);

private:
    FPType r_;
    FPType i_;
};

template<typename FPType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator=(const FPType& __t) {
    r_ = __t;
    i_ = FPType();
    return *this;
}

template<typename FPType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator*=(const FPType& __t) {
    r_ *= __t;
    i_ *= __t;
    return *this;
}

template<typename FPType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator/=(const FPType& __t) {
    r_ /= __t;
    i_ /= __t;
    return *this;
}

template<typename FPType>
template<typename CXType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator=(const complex<CXType>& __z) {
    r_ = __z.real();
    i_ = __z.imag();
    return *this;
}

template<typename FPType>
template<typename CXType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator+=(const complex<CXType>& __z) {
    r_ += __z.real();
    i_ += __z.imag();
    return *this;
}

template<typename FPType>
template<typename CXType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator-=(const complex<CXType>& __z) {
    r_ -= __z.real();
    i_ -= __z.imag();
    return *this;
}

template<typename FPType>
template<typename CXType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator*=(const complex<CXType>& __z) {
    const FPType __r = r_ * __z.real() - i_ * __z.imag();
    i_ = r_ * __z.imag() + i_ * __z.real();
    r_ = __r;
    return *this;
}

template<typename FPType>
template<typename CXType>
SYCL_EXTERNAL complex<FPType>& complex<FPType>::operator/=(const complex<CXType>& __z) {
    const FPType __r =  r_ * __z.real() + i_ * __z.imag();
    const FPType __x = __z.real();
    const FPType __y = __z.imag();
    const FPType __n = __x * __x + __y * __y;
    i_ = (i_ * __z.real() - r_ * __z.imag()) / __n;
    r_ = __r / __n;
    return *this;
}

// Operators:
//@{
///  Return new complex value @a x plus @a y.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator+(const complex<FPType>& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __x;
    __r += __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator+(const complex<FPType>& __x, const FPType& __y) {
    complex<FPType> __r = __x;
    __r += __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator+(const FPType& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __y;
    __r += __x;
    return __r;
}
//@}
//@{
///  Return new complex value @a x minus @a y.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator-(const complex<FPType>& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __x;
    __r -= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator-(const complex<FPType>& __x, const FPType& __y) {
    complex<FPType> __r = __x;
    __r -= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator-(const FPType& __x, const complex<FPType>& __y) {
    complex<FPType> __r(__x, -__y.imag());
    __r -= __y.real();
    return __r;
}
//@}
//@{
///  Return new complex value @a x times @a y.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator*(const complex<FPType>& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __x;
    __r *= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator*(const complex<FPType>& __x, const FPType& __y) {
    complex<FPType> __r = __x;
    __r *= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator*(const FPType& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __y;
    __r *= __x;
    return __r;
}
//@}
//@{
///  Return new complex value @a x divided by @a y.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator/(const complex<FPType>& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __x;
    __r /= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator/(const complex<FPType>& __x, const FPType& __y) {
    complex<FPType> __r = __x;
    __r /= __y;
    return __r;
}

template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator/(const FPType& __x, const complex<FPType>& __y) {
    complex<FPType> __r = __x;
    __r /= __y;
    return __r;
}
//@}
///  Return @a x.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator+(const complex<FPType>& __x) { return __x;
}

///  Return complex negation of @a x.
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> operator-(const complex<FPType>& __x) {
    return complex<FPType>(-__x.real(), -__x.imag());
}

// Return complex conjugate
template<typename FPType>
SYCL_EXTERNAL inline complex<FPType> conj(const complex<FPType>& __z) {
    return complex<FPType>(__z.real(), -__z.imag());
}

}
#endif
