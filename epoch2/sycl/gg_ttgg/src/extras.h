#ifndef EXTRAS_H
#define EXTRAS_H 1

#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif

namespace extras {

template< typename FPType >
class complex {
public:

  SYCL_EXTERNAL complex(const complex<FPType> &c) : r_(c.r_), i_(c.i_) { }
  SYCL_EXTERNAL complex<FPType>& operator=(const complex<FPType> &c) {
    if (this != &c) {
      r_ = c.r_;
      i_ = c.i_;
    }
    return *this;
  }
  SYCL_EXTERNAL ~complex() { }

  SYCL_EXTERNAL complex() : r_(0.), i_(0.) { }
  SYCL_EXTERNAL complex(FPType r, FPType i) : r_(r), i_(i) { }
  SYCL_EXTERNAL complex(FPType r) : r_(r), i_(0) { }

  SYCL_EXTERNAL FPType real() const { return r_; }
  SYCL_EXTERNAL FPType imag() const { return i_; }

  SYCL_EXTERNAL complex<FPType>& operator/(const complex<FPType> &z2) {
    const FPType s1 = (this->real() * z2.real())+(this->imag() * z2.imag());
    const FPType s2 = (this->imag()*z2.real()) - (this->real() * z2.imag());
    const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
    r_ = s1/s3;
    i_ = s2/s3;
    return *this;
  }

  SYCL_EXTERNAL const complex<FPType> operator/(const complex<FPType> &z2) const {
    const FPType s1 = (this->real() * z2.real())+(this->imag() * z2.imag());
    const FPType s2 = (this->imag()*z2.real()) - (this->real() * z2.imag());
    const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
    return complex<FPType>(s1/s3,s2/s3);
  }

  SYCL_EXTERNAL complex<FPType>& operator*(const complex<FPType> &z2) {
    const FPType s1 = this->real() * z2.real();
    const FPType s2 = this->imag() * z2.imag();
    const FPType s3 = (this->real() + this->imag())*(z2.real() + z2.imag());
    r_ = s1-s2;
    i_ = s3-s1-s2;
    return *this;
  }

  SYCL_EXTERNAL const complex<FPType> operator*(const complex<FPType> &z2) const {
    const FPType s1 = this->real() * z2.real();
    const FPType s2 = this->imag() * z2.imag();
    const FPType s3 = (this->real() + this->imag())*(z2.real() + z2.imag());
    return complex<FPType>(s1-s2,s3-s1-s2);
  }

  SYCL_EXTERNAL complex<FPType>& operator+(const complex<FPType> &z2) {
    r_ += z2.real();
    i_ += z2.imag();
    return *this;
  }

  SYCL_EXTERNAL const complex<FPType> operator+(const complex<FPType> &z2) const {
    return complex<FPType>(*this) + z2;
  }

  SYCL_EXTERNAL complex<FPType>& operator+=(const complex<FPType> &z2) {
    r_ += z2.real();
    i_ += z2.imag();
    return *this;
  }

  SYCL_EXTERNAL complex<FPType>& operator-(const complex<FPType> &z2) {
    r_ -= z2.real();
    i_ -= z2.imag();
    return *this;
  }

  SYCL_EXTERNAL const complex<FPType> operator-(const complex<FPType> &z2) const {
    return complex<FPType>(*this) - z2;
  }

  FPType r_;
  FPType i_;
};

}
#endif
