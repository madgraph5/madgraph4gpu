#include <cfenv> // for fetestexcept
#include <cfloat> // for FLT_MIN
#include <cmath> // for sqrt
#include <iostream>

void printFPEconfig()
{
  int fpes = fegetexcept();
  std::cout << "Floating Point Exception: analyse fegetexcept()=" << fpes << std::endl;
  std::cout << "  FPE trap for FE_DIVBYZERO is" << ( ( fpes & FE_DIVBYZERO ) ? " " : " NOT " ) << "enabled" << std::endl;
  std::cout << "  FPE trap for FE_INEXACT is" << ( ( fpes & FE_INEXACT ) ? " " : " NOT " ) << "enabled" << std::endl;
  std::cout << "  FPE trap for FE_INVALID is" << ( ( fpes & FE_INVALID ) ? " " : " NOT " ) << "enabled" << std::endl;
  std::cout << "  FPE trap for FE_OVERFLOW is" << ( ( fpes & FE_OVERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
  std::cout << "  FPE trap for FE_UNDERFLOW is" << ( ( fpes & FE_UNDERFLOW ) ? " " : " NOT " ) << "enabled" << std::endl;
}

void printFPEs()
{
  if( std::fetestexcept( FE_ALL_EXCEPT ) )
    std::cerr << "Floating Point Exception: FE_ALL_EXCEPT reported" << std::endl;
  else
    std::cerr << "Floating Point Exception: no FPEs reported (FE_ALL_EXCEPT not reported)" << std::endl;
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
  //if( std::fetestexcept( FE_DENORMAL ) ) // this does not exist...
  //  std::cerr << "Floating Point Exception: FE_DENORMAL reported" << std::endl;
}

void clearFPEs()
{
  std::cerr << "Floating Point Exception: clear all exceptions" << std::endl;
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
  void testfpes_( const bool* pClearFPEs )
  {
    printFPEconfig();

    std::cout << std::endl;
    printFPEs();

    std::cout << std::endl;
    testUnderflowF();
    if( *pClearFPEs ) clearFPEs();

    std::cout << std::endl;
    testUnderflowD();
    if( *pClearFPEs ) clearFPEs();

    std::cout << std::endl;
    testOverflow();
    if( *pClearFPEs ) clearFPEs();

    std::cout << std::endl;
    testInvalid();
    if( *pClearFPEs ) clearFPEs();

    std::cout << std::endl;
    testDivBy0();
    if( *pClearFPEs ) clearFPEs();

    std::cout << std::endl;
    printFPEs();
 }
}
