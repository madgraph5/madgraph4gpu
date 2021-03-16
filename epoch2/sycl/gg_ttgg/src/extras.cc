#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif
#include <math.h>
#include "extras.h"

namespace extras
{
  
double complex::real(void) const{
  return r;
}

double complex::imag(void) const{
  return i;
}

double complex::abs(void){
  return sycl::sqrt(r*r+i*i);
}

double complex::arg(void){
  double phi;
  if(i==0)
    phi = 90;
  else
    phi = sycl::atan(i/r)*180/M_PI;
  return phi;
}

void complex::DivError(void){
  div_zero = true;
}

bool complex::DivByZero(void){
  bool flag = div_zero;
  div_zero = false;
  return flag;
}

void complex::setpolar(double mag, double phi)
{
  phi = phi * M_PI / 180;
  r = mag*sycl::cos(phi);
  i = mag*sycl::sin(phi);
}

complex operator+(complex z1,complex z2){
  double r_temp = z1.real() + z2.real();
  double i_temp= z1.imag() + z2.imag();
  complex temp(r_temp, i_temp);
  return temp;
}

complex operator+=(complex z1,complex z2){
  return operator+(z1, z2);
}

complex operator+(complex z2){
  complex z1(1.0);
  return z1 * z2;
}

complex operator-(complex z2){
  complex z1(-1.0);
  return z1 * z2;
}

complex operator-(complex z1, complex z2){
  double r_temp = z1.real() - z2.real();
  double i_temp = z1.imag() - z2.imag();
  complex temp(r_temp, i_temp);
  return temp;
}

complex operator/(complex z1,complex z2){
  complex temp;

  double mag2 = z2.abs();
  if(mag2==0.0){
    temp.DivError();
  }
  else{
    double mag1 = z1.abs();
    double phi = z1.arg();
    mag1 /= mag2;
    phi -= z2.arg();
    temp.setpolar(mag1, phi);
  }
  return temp;
}

complex operator*(complex z1,complex z2){
  complex temp;
  double mag = z1.abs();
  mag *= z2.abs();
  double phi = z1.arg();
  phi += z2.arg();
  temp.setpolar(mag,phi);
  return temp;
}

complex operator*(double r ,complex z2){
  complex z1 = complex(r);
  return z1 * z2;
}

complex operator/(complex z1, double r){
  complex z2=complex(r);
  complex temp = z1/z2;
  return temp;
}

complex conj(const complex& z){
  return complex(z.real(), -z.imag()); 
}

}  //end of namespace
