#include <ompnumthreads.h>

// NB1: The C function ompnumthreadsNotSetMeansOneThread_ is called by Fortran code
// Hence the trailing "_": 'call xxx()' links to xxx_
// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

// NB2: This file also contains C++ code and is built using g++
// Hence use 'extern "C"' to avoid name mangling by the C++ compiler
// See https://www.geeksforgeeks.org/extern-c-in-c

extern "C"
{
  void ompnumthreads_not_set_means_one_thread_()
  {
    const int debuglevel = 0;                        // quiet(-1), info(0), debug(1)
    ompnumthreadsNotSetMeansOneThread( debuglevel ); // call the inline C++ function defined in the .h file
  }
}
