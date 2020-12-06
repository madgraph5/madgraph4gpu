#ifndef ALSIMPLE_H
#define ALSIMPLE_H 1

#include <random>

namespace alsimple {

//--------------------------------------------------------------------------

enum alsrandType_t {
  ALSIMPLE_RNG_PSEUDO_MT19937
};

enum randStatus_t {
  RANDOK = 0,
  RANDBAD
};

class randGenerator {
  public:

  randGenerator();
  ~randGenerator() { }

  std::mt19937 gen;
  std::uniform_real_distribution<double> dis_double;
  std::uniform_real_distribution<float> dis_float;
};
typedef randGenerator* randGenerator_t;

randStatus_t createGeneratorHost( randGenerator_t *gen, alsrandType_t type );

randStatus_t createGenerator( randGenerator_t *gen, alsrandType_t type );

randStatus_t randSetPseudoRandomGeneratorSeed( randGenerator_t gen, unsigned long long seed );

randStatus_t randDestroyGenerator( randGenerator_t gen );

randStatus_t randGenerateUniformDouble( randGenerator_t gen, double *out, size_t nsize);

randStatus_t randGenerateUniform( randGenerator_t gen, float *out, size_t nsize);

//--------------------------------------------------------------------------


template< typename FPType >
class complex {
public:

  complex() : r_(0.), i_(0.) { }
  complex(FPType r, FPType i) : r_(r), i_(i) { }
  ~complex() { }

  FPType real() const { return r_; }
  FPType imag() const { return i_; }

  complex<FPType>& operator/(const complex<FPType> &z2) {
    const FPType s1 = (this->real() * z2.real())+(this->imag() * z2.imag());
    const FPType s2 = (this->imag()*z2.real()) - (this->real() * z2.imag());
    const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
    r_ = s1/s3;
    i_ = s2/s3;
    return *this;
  }

  const complex<FPType> operator/(const complex<FPType> &z2) const {
    return complex<FPType>(*this) / z2;
  }

  complex<FPType>& operator*(const complex<FPType> &z2) {
    const FPType s1 = this->real() * z2.real();
    const FPType s2 = this->imag() * z2.imag();
    const FPType s3 = (this->real() + this->imag())*(z2.real() + z2.imag());
    r_ = s1-s2;
    i_ = s3-s1-s2;
    return *this;
  }

  const complex<FPType> operator*(const complex<FPType> &z2) const {
    return complex<FPType>(*this) * z2;
  }

  complex<FPType>& operator+(const complex<FPType> &z2) {
    r_ += z2.real();
    i_ += z2.imag();
    return *this;
  }

  const complex<FPType> operator+(const complex<FPType> &z2) const {
    return complex<FPType>(*this) + z2;
  }

  complex<FPType>& operator-(const complex<FPType> &z2) {
    r_ -= z2.real();
    i_ -= z2.imag();
    return *this;
  }

  const complex<FPType> operator-(const complex<FPType> &z2) const {
    return complex<FPType>(*this) - z2;
  }

  FPType r_;
  FPType i_;
};

}
#endif
