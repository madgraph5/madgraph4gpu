#pragma once
#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif

namespace extras
{
class complex {
  public:
    complex() :	r(0.0), i(0.0) {}

    complex(const complex& z) : r(z.real()), i(z.imag()) {} // copy constructor

    complex(double a, double b) : r(a), i(b){}

    complex(double a) : r(a), i(0.0){}
    void SetPolar(double, double); // second argument in degrees
    SYCL_EXTERNAL
    double real(void) const;
    SYCL_EXTERNAL
    double imag(void) const;
    double abs(void);
    double arg(void); // in degrees
    void setpolar(double, double); // second argument in degree
    void DivError(void); // used by divide operator
    bool DivByZero(void);
  private:
    double r;
    double i;
    bool div_zero;
};
SYCL_EXTERNAL
complex operator+(complex);
SYCL_EXTERNAL
complex operator+=(complex, complex);
SYCL_EXTERNAL
complex operator-(complex);
SYCL_EXTERNAL
complex operator+(complex,complex);
SYCL_EXTERNAL
complex operator-(complex,complex);
SYCL_EXTERNAL
complex operator*(complex,complex);
SYCL_EXTERNAL
complex operator*(double ,complex);
SYCL_EXTERNAL
complex operator/(complex,complex);
SYCL_EXTERNAL
complex operator/(complex, double);
SYCL_EXTERNAL
complex conj(const complex& z); 

}  //end of namespace
