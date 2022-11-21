#ifndef ALSIMPLE_H
#define ALSIMPLE_H 1

#include "alpaka/alpaka.hpp"

namespace alsimple {


template< typename FPType >
class complex {
public:

  ALPAKA_FN_ACC complex(const complex<FPType> &c) : r_(c.r_), i_(c.i_) { }
  ALPAKA_FN_ACC complex<FPType>& operator=(const complex<FPType> &c) {
    if (this != &c) {
      r_ = c.r_;
      i_ = c.i_;
    }
    return *this;
  }
  //ALPAKA_FN_ACC ~complex() { }

  ALPAKA_FN_ACC complex() { }
  ALPAKA_FN_ACC complex(FPType r, FPType i) : r_(r), i_(i) { }

  ALPAKA_FN_ACC FPType real() const { return r_; }
  ALPAKA_FN_ACC FPType imag() const { return i_; }

  ALPAKA_FN_ACC complex<FPType>& operator/(const complex<FPType> &z2) {
    const FPType s1 = (this->real() * z2.real())+(this->imag() * z2.imag());
    const FPType s2 = (this->imag()*z2.real()) - (this->real() * z2.imag());
    const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
    r_ = s1/s3;
    i_ = s2/s3;
    return *this;
  }

  ALPAKA_FN_ACC const complex<FPType> operator/(const complex<FPType> &z2) const {
//    return complex<FPType>(*this) / z2;
    const FPType s1 = (this->real() * z2.real())+(this->imag() * z2.imag());
    const FPType s2 = (this->imag()*z2.real()) - (this->real() * z2.imag());
    const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
    return complex<FPType>(s1/s3,s2/s3);
  }

  ALPAKA_FN_ACC complex<FPType>& operator*(const complex<FPType> &z2) {
    const FPType s1 = this->real() * z2.real();
    const FPType s2 = this->imag() * z2.imag();
    const FPType s3 = (this->real() + this->imag())*(z2.real() + z2.imag());
    r_ = s1-s2;
    i_ = s3-s1-s2;
    return *this;
  }

  ALPAKA_FN_ACC const complex<FPType> operator*(const complex<FPType> &z2) const {
//    return complex<FPType>(*this) * z2;
    const FPType s1 = this->real() * z2.real();
    const FPType s2 = this->imag() * z2.imag();
    const FPType s3 = (this->real() + this->imag())*(z2.real() + z2.imag());
    return complex<FPType>(s1-s2,s3-s1-s2);
  }

  ALPAKA_FN_ACC complex<FPType>& operator+(const complex<FPType> &z2) {
    r_ += z2.real();
    i_ += z2.imag();
    return *this;
  }

  ALPAKA_FN_ACC const complex<FPType> operator+(const complex<FPType> &z2) const {
    return complex<FPType>(*this) + z2;
  }

  ALPAKA_FN_ACC complex<FPType>& operator-(const complex<FPType> &z2) {
    r_ -= z2.real();
    i_ -= z2.imag();
    return *this;
  }

  ALPAKA_FN_ACC const complex<FPType> operator-(const complex<FPType> &z2) const {
    return complex<FPType>(*this) - z2;
  }

  ALPAKA_FN_ACC complex<FPType>& operator-=(const complex<FPType> &z2) {
    r_ -= z2.real();
    i_ -= z2.imag();
    return *this;
  }

  ALPAKA_FN_ACC complex<FPType>& operator+=(const complex<FPType> &z2) {
    r_ += z2.real();
    i_ += z2.imag();
    return *this;
  }

  FPType r_;
  FPType i_;
};

}
#endif
