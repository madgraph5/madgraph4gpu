#include <cfenv> // for fetestexcept
#include <cfloat> // for FLT_MIN
#include <cmath> // for sqrt
#include <iostream>

void printFPEs()
{
  if( std::fetestexcept( FE_ALL_EXCEPT ) )
    std::cerr << "Floating Point Exception: FE_ALL_EXCEPT reported" << std::endl;
  else
    std::cerr << "Floating Point Exception: FE_ALL_EXCEPT NOT reported" << std::endl;
  if( std::fetestexcept( FE_DIVBYZERO ) )
    std::cerr << "Floating Point Exception: FE_DIVBYZERO reported" << std::endl;
  if( std::fetestexcept( FE_INVALID ) )
    std::cerr << "Floating Point Exception: FE_INVALID reported" << std::endl;
  if( std::fetestexcept( FE_OVERFLOW ) )
    std::cerr << "Floating Point Exception: FE_OVERFLOW reported" << std::endl;
  if( std::fetestexcept( FE_UNDERFLOW ) )
    std::cerr << "Floating Point Exception: FE_UNDERFLOW reported" << std::endl;
  if( std::fetestexcept( FE_INEXACT ) ) // this should not throw a signal?
    std::cerr << "Floating Point Exception: FE_INEXACT reported" << std::endl;
  std::feclearexcept( FE_ALL_EXCEPT );  
}

void testUnderflowF()
{
  std::cout << "UNDERFLOW" << std::endl;
  std::cout << "FLT_MIN = " << FLT_MIN << std::endl;
  float underf1 = FLT_MIN;
  underf1 /= (float)10;
  std::cout << "FLT_MIN / 10 = " << underf1 << std::endl;
  float underf2 = FLT_MIN;
  underf2 *= underf2;
  std::cout << "FLT_MIN**2 = " << underf2 << std::endl;
  printFPEs();
}

void testUnderflowD()
{
  std::cout << "UNDERFLOW" << std::endl;
  std::cout << "DBL_MIN = " << DBL_MIN << std::endl;
  double underd1 = DBL_MIN;
  underd1 /= (double)10;
  std::cout << "DBL_MIN / 10 = " << underd1 << std::endl;
  double underd2 = DBL_MIN;
  underd2 *= underd2;
  std::cout << "DBL_MIN**2 = " << underd2 << std::endl;
  printFPEs();
}

void testOverflow()
{
  std::cout << "OVERFLOW" << std::endl;
  std::cout << "FLT_MAX = " << FLT_MAX << std::endl;
  float overf = FLT_MAX * (float)10;
  std::cout << "FLT_MAX * 10 = " << overf << std::endl;
  std::cout << "DBL_MAX = " << DBL_MAX << std::endl;
  double overd = DBL_MAX * (double)10;
  std::cout << "DBL_MAX * 10 = " << overd << std::endl;
  printFPEs();
}

void testInvalid()
{
  std::cout << "INVALID" << std::endl;
  double cI = sqrt( -1. );
  std::cout << "sqrt( -1. ) = " << cI << std::endl;
  printFPEs();
}

void testDivBy0()
{
  std::cout << "DIV_BY_ZERO" << std::endl;
  double oneBy0 = 1. / 0.;
  std::cout << "1. / 0. = " << oneBy0 << std::endl;
  printFPEs();
}

extern "C"
{
  void testfpes_()
  {
    printFPEs();
    testUnderflowF();
    std::cout << std::endl;
    testUnderflowD();
    std::cout << std::endl;
    testOverflow();
    std::cout << std::endl;
    testInvalid();
    std::cout << std::endl;
    testDivBy0();
  }
}
