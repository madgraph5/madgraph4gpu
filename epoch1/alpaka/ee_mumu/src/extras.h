#pragma once

#include <cuda_to_cupla.hpp>

#include <iostream>

namespace extras
{

class extrasGenerator_t {
public:
  extrasGenerator_t() {
    setseed(0x249d070a);
  }
  void uniform(double *dout, const size_t n) {
    /*     -----------------
     * universal random number generator proposed by marsaglia and zaman
     * in report fsu-scri-87-50
     * in this version rvec is a double precision variable. */
    for(size_t itr = 0; itr<n; ++itr) {
      double uni = ranu[iranmr] - ranu[jranmr];
      if (uni < 0)
        uni = uni + 1;
      ranu[iranmr] = uni;
      iranmr = iranmr - 1;
      jranmr = jranmr - 1;
      if (iranmr == 0)
        iranmr = 97;
      if (jranmr == 0)
        jranmr = 97;
      ranc = ranc - rancd;
      if (ranc < 0)
        ranc = ranc + rancm;
      uni = uni - ranc;
      if (uni < 0)
        uni = uni + 1;
      if (uni<= 1e-16) {
        // redo
        itr--;
        continue;
      }
//dhs
//std::cout << "dout[" << itr << "]=" << uni << std::endl;
      dout[itr] = uni;
    }
  }

  void setseed(unsigned long long seed) {
    uint32_t x = seed ^ (seed>>32);
    int ij = (x&65535)%31329;
    int kl = (x>>16)%30082;
//dhs
//std::cout << "setseed ij=" << ij << " kl=" << kl << std::endl;
    /*     -----------------
     * initializing routine for ranmar, must be called before generating
     * any pseudorandom numbers with ranmar. the input values should be in
     * the ranges 0<=ij<=31328 ; 0<=kl<=30081 */
    /* this shows correspondence between the simplified input seeds ij, kl
     * and the original marsaglia-zaman seeds i,j,k,l.
     * to get the standard values in the marsaglia-zaman paper (i=12,j=34
     * k=56,l=78) put ij=1802, kl=9373 */
    int i = ij / 177 % 177 + 2;
    int j = ij % 177 + 2;
    int k = (kl / 169) % 178 + 1;
    int l = kl % 169;
    for (int ii = 1; ii < 98; ii++) {
      double s = 0;
      double t = .5;
      for (int jj = 1; jj < 25; jj++) {
        int m = ((i * j % 179) * k) % 179;
        i = j;
        j = k;
        k = m;
        l = (53 * l + 1) % 169;
        if ((l * m) % 64 >= 32)
          s = s + t;
        t = .5 * t;
      }
      ranu[ii] = s;
    }
    ranc = 362436. / 16777216.;
    rancd = 7654321. / 16777216.;
    rancm = 16777213. / 16777216.;
    iranmr = 97;
    jranmr = 33;
  }
  static void rnCreate(extrasGenerator_t **p) { *p = new extrasGenerator_t; return; }
  static void rnDestroy(extrasGenerator_t *p) { delete p; return; }

private:
  double ranu[98];
  double ranc, rancd, rancm;
  int iranmr, jranmr;
};

// assumed to be a pod of exactly 2* sizeof(fptype)
template<typename FPType>
class cxtypeparam {
  public:
    FPType r;
    FPType i;
};

// forward decl
template<typename T_Acc, typename FPType>
class cxtypeop;

template<typename T_Acc, typename FPType>
class cxtyperef {
  public:
    cxtyperef() = delete;
    cxtyperef(const T_Acc &acc, cxtypeparam<FPType> &cr) : acc_(acc), cp_(cr) { }
    cxtyperef &operator=(const cxtypeparam<FPType> &other) {
      this->cp_ = other;
      return *this;
    }
    cxtyperef &operator=(const cxtyperef &other) {
      this->cp_ = other.cp_;
      return *this;
    }
    cxtyperef &operator=(const cxtypeop<T_Acc, FPType> &other) {
      this->cp_ = other.c_;
      return *this;
    }

    FPType real(void) const { return cp_.r; }
    FPType imag(void) const { return cp_.i; }

    const T_Acc &acc_;
    cxtypeparam<FPType> &cp_;
};

template<typename T_Acc, typename FPType>
class cxtypeconstref {
  public:
    cxtypeconstref() = delete;
    cxtypeconstref(const T_Acc &acc, const cxtypeparam<FPType> &cr) : acc_(acc), cp_(cr) { }

    FPType real(void) const { return cp_.r; }
    FPType imag(void) const { return cp_.i; }

    const T_Acc &acc_;
    const cxtypeparam<FPType> &cp_;
};

template<typename T_Acc, typename FPType>
class cxtypeptr {
  public:
    cxtypeptr() = delete;
    cxtypeptr(const T_Acc &acc, cxtypeparam<FPType> *cr) : acc_(acc), cp_(cr) { }
    cxtypeptr(const T_Acc &acc, const cxtypeparam<FPType> *cr) : acc_(acc), cp_(cr) { }
    cxtypeptr &operator=(const cxtypeparam<FPType> *other) {
      this->cp_ = other;
      return *this;
    }
    cxtypeptr &operator=(const cxtypeptr &other) {
      this->cp_ = other.cp_;
      return *this;
    }
    cxtypeptr &operator=(const cxtypeop<T_Acc,FPType> &other) {
      *this->cp_ = other.c_;
      return *this;
    }
    cxtyperef<T_Acc,FPType> &operator*() {
      return cxtyperef<T_Acc,FPType>(acc_, *cp_);
    }
    cxtyperef<T_Acc,FPType> operator[](const uint i) {
      return cxtyperef<T_Acc,FPType>(acc_, cp_[i]);
    }

    const T_Acc &acc_;
    cxtypeparam<FPType> *cp_;
};

template<typename T_Acc, typename FPType>
class cxtypeop {
  public:
    cxtypeop() = delete;
    cxtypeop(const T_Acc &acc, const cxtypeparam<FPType> &cr) : acc_(acc), c_(cr) { };
    cxtypeop(const cxtyperef<T_Acc,FPType> &cr) : acc_(cr.acc_), c_(cr.cp_) { }
    cxtypeop(const cxtypeconstref<T_Acc,FPType> &cr) : acc_(cr.acc_), c_(cr.cp_) { }
    cxtypeop(const T_Acc &acc) : acc_(acc) { c_ = { .r = 0., .i = 0. }; };
    cxtypeop(const T_Acc &acc, const FPType r) : acc_(acc) { c_ = { .r = r, .i = 0. }; };
    cxtypeop(const T_Acc &acc, const FPType r, const FPType i) : acc_(acc) { c_ = { .r= r, .i=i }; };

    FPType real(void) const;
    FPType imag(void) const;
    FPType abs(void) const;
    FPType arg(void) const; // in degrees
    void setpolar(FPType, FPType); // second argument in degree
    void DivError(void); // used by divide operator
    bool DivByZero(void);

    const T_Acc &acc_;
    cxtypeparam<FPType> c_;
};

template<typename T_Acc, typename FPType> FPType cxtypeop<T_Acc,FPType>::real(void) const{
  return c_.r;
}

template<typename T_Acc, typename FPType> FPType cxtypeop<T_Acc,FPType>::imag(void) const{
  return c_.i;
}

template<typename T_Acc, typename FPType> FPType cxtypeop<T_Acc,FPType>::abs(void) const{
  return alpaka::math::sqrt(acc_, c_.r*c_.r+c_.i*c_.i);
}

template<typename T_Acc, typename FPType> FPType cxtypeop<T_Acc,FPType>::arg(void) const{
  FPType phi;
  if(c_.r==0) {
    if (c_.i>0) {
      phi = 90.;
    } else {
      phi = -90.;
    }
  } else {
    if (c_.r > 0) {
      phi = alpaka::math::atan(acc_, c_.i/c_.r)*180./M_PI;
     } else {
      if (c_.i >=0) {
        phi = alpaka::math::atan(acc_, c_.i/c_.r)*180./M_PI + 180.;
      } else {
        phi = alpaka::math::atan(acc_, c_.i/c_.r)*180./M_PI - 180.;
      }
     }
  }
  return phi;
}

template<typename T_Acc, typename FPType> void cxtypeop<T_Acc,FPType>::DivError(void){
//  c_.div_zero = true;
}

template<typename T_Acc, typename FPType> bool cxtypeop<T_Acc,FPType>::DivByZero(void){
  return false;
//  bool flag = c_.div_zero;
//  c_.div_zero = false;
//  return flag;
}

template<typename T_Acc, typename FPType> void cxtypeop<T_Acc,FPType>::setpolar(FPType mag, FPType phi)
{
  phi = phi * M_PI / 180;
  c_.r = mag*alpaka::math::cos(acc_, phi);
  c_.i = mag*alpaka::math::sin(acc_, phi);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const FPType r, const cxtypeop<T_Acc,FPType> &z2){
  cxtypeop<T_Acc,FPType> temp(z2.acc_,r);
  return temp + z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const FPType r, const cxtyperef<T_Acc,FPType> &z2){
  return r + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const FPType r, const cxtypeconstref<T_Acc,FPType> &z2){
  return r + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeop<T_Acc,FPType> &z1, const cxtyperef<T_Acc,FPType> &z2) {
  return z1 + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeop<T_Acc,FPType> &z1, const cxtypeconstref<T_Acc,FPType> &z2) {
  return z1 + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtyperef<T_Acc,FPType> &z1, const cxtypeop<T_Acc,FPType> &z2) {
  return cxtypeop<T_Acc,FPType>(z1) + z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeconstref<T_Acc,FPType> &z1, const cxtypeop<T_Acc,FPType> &z2) {
  return cxtypeop<T_Acc,FPType>(z1) + z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtyperef<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) + cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeop<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  FPType r_temp = z1.real() + z2.real();
  FPType i_temp= z1.imag() + z2.imag();
  cxtypeop<T_Acc,FPType> temp(z1.acc_, r_temp, i_temp);
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeop<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return z1 - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeop<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return z1 - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtyperef<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) - z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) - z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const FPType r, const cxtypeop<T_Acc,FPType> &z2){
  cxtypeop<T_Acc,FPType> temp(z2.acc_,r);
  return temp - z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const FPType r, const cxtyperef<T_Acc,FPType> &z2){
  return r - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const FPType r, const cxtypeconstref<T_Acc,FPType> &z2){
  return r - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtyperef<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) - cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeop<T_Acc,FPType> &z1, const cxtypeop<T_Acc,FPType> &z2){
  FPType r_temp = z1.real() - z2.real();
  FPType i_temp = z1.imag() - z2.imag();
  cxtypeop<T_Acc,FPType> temp(z1.acc_,r_temp, i_temp);
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeop<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return z1 / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeop<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return z1 / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtyperef<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) / z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) / z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtyperef<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeop<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
/*
  cxtypeop<T_Acc,FPType> temp(z1);

  FPType mag2 = z2.abs();
  if(mag2==0.0){
//    temp.DivError();
  }
  else{
    FPType mag1 = z1.abs();
    FPType phi = z1.arg();
    mag1 /= mag2;
    phi -= z2.arg();
    temp.setpolar(mag1, phi);
  }
  return temp;
*/
  const FPType s1 = (z1.real() * z2.real())+(z1.imag() * z2.imag());
  const FPType s2 = (z1.imag()*z2.real()) - (z1.real() * z2.imag());
  const FPType s3 = (z2.real()*z2.real()) + (z2.imag() * z2.imag());
  cxtypeop<T_Acc,FPType> temp(z1.acc_);
  temp.c_.r = s1 / s3;
  temp.c_.i = s2 / s3;
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeop<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return z1 * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeop<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return z1 * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtyperef<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtyperef<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeconstref<T_Acc,FPType> &z1,const cxtyperef<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtyperef<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeconstref<T_Acc,FPType> &z1,const cxtypeconstref<T_Acc,FPType> &z2){
  return cxtypeop<T_Acc,FPType>(z1) * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeop<T_Acc,FPType> &z1,const cxtypeop<T_Acc,FPType> &z2){
//  cxtypeop<T_Acc,FPType> temp(z1);
//  FPType mag = z1.abs();
//  mag *= z2.abs();
//  FPType phi = z1.arg();
//  phi += z2.arg();
//  temp.setpolar(mag,phi);
//  return temp;
  const FPType s1 = z1.real() * z2.real();
  const FPType s2 = z1.imag() * z2.imag();
  const FPType s3 = (z1.real() + z1.imag())*(z2.real() + z2.imag());
  cxtypeop<T_Acc,FPType> temp(z1.acc_);
  temp.c_.r = s1-s2;
  temp.c_.i = s3-s1-s2;
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const FPType r,const cxtyperef<T_Acc,FPType> &z2){
  return r * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const FPType r,const cxtypeconstref<T_Acc,FPType> &z2){
  return r * cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const FPType r ,const cxtypeop<T_Acc,FPType> &z2){
  cxtypeop<T_Acc,FPType> z1(z2.acc_,r);
  return z1 * z2;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtyperef<T_Acc,FPType> &z1, const FPType r){
  return r * cxtypeop<T_Acc,FPType>(z1);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeconstref<T_Acc,FPType> &z1, const FPType r){
  return r * cxtypeop<T_Acc,FPType>(z1);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator*(const cxtypeop<T_Acc,FPType> &z1, const FPType r){
  return r * z1;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtyperef<T_Acc,FPType> &z){
  return -cxtypeop<T_Acc,FPType>(z);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeconstref<T_Acc,FPType> &z){
  return -cxtypeop<T_Acc,FPType>(z);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator-(const cxtypeop<T_Acc,FPType> &z){
  return -1.0*z;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtyperef<T_Acc,FPType> &z){
  return +cxtypeop<T_Acc,FPType>(z);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeconstref<T_Acc,FPType> &z){
  return +cxtypeop<T_Acc,FPType>(z);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator+(const cxtypeop<T_Acc,FPType> &z){
  return z;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtyperef<T_Acc,FPType> &z1,const FPType r){
  return cxtypeop<T_Acc,FPType>(z1) / r;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeconstref<T_Acc,FPType> &z1,const FPType r){
  return cxtypeop<T_Acc,FPType>(z1) / r;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const cxtypeop<T_Acc,FPType> &z1, const FPType r){
  cxtypeop<T_Acc,FPType> z2(z1.acc_,r);
  cxtypeop<T_Acc,FPType> temp = z1/z2;
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const FPType r, const cxtyperef<T_Acc,FPType> &z2){
  return r / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const FPType r, const cxtypeconstref<T_Acc,FPType> &z2){
  return r / cxtypeop<T_Acc,FPType>(z2);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> operator/(const FPType r, const cxtypeop<T_Acc,FPType> &z2) {
  cxtypeop<T_Acc,FPType> z1(z2.acc_,r);
  cxtypeop<T_Acc,FPType> temp = z1/z2;
  return temp;
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> conj(const cxtyperef<T_Acc,FPType> &z) {
  return cxtypeop<T_Acc,FPType>(z.acc_, z.cp_.r, -z.cp_.i);
}

template<typename T_Acc, typename FPType> cxtypeop<T_Acc,FPType> conj(const cxtypeconstref<T_Acc,FPType> &z) {
  return cxtypeop<T_Acc,FPType>(z.acc_, z.cp_.r, -z.cp_.i);
}

}  //end of namespace
