#pragma once

namespace extras
{
class complex {
  public:
    complex() :	r(0.0), i(0.0) {}
    complex(const complex& z) : r(z.real()), i(z.imag()) {} // copy constructor
    complex(double a, double b) : r(a), i(b){}
    complex(double a) : r(a), i(0.0){}
    void SetPolar(double, double); // second argument in degrees
    double real(void) const;
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

complex operator+(complex);
complex operator-(complex);
complex operator+(complex,complex);
complex operator-(complex,complex);
complex operator*(complex,complex);
complex operator*(double ,complex);
complex operator/(complex,complex);
complex operator/(complex, double);

complex conj(const complex& z); 

}  //end of namespace
