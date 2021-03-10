#include "simpletimermap.h"

#include <fstream>
#include <stdio.h>

// NB1: The C functions counters_xxx_ in this file are called by Fortran code
// Hence the trailing "_": 'call counters_end()' links to counters_end_
// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

// NB2: This file also contains C++ code and is built using g++
// Hence use 'extern "C"' to avoid name mangling by the C++ compiler
// See https://www.geeksforgeeks.org/extern-c-in-c

extern "C"
{

  static mgOnGpu::SimpleTimerMap timermap;

  static int counters_counter1 = 0;
  static int counters_counter2 = 0;

  void counters_initialise_()
  {
    return;
  }

  void counters_start_( char* key )
  {
    counters_counter1++;
    timermap.start( std::string( key ) );
    return;
  }

  void counters_end_()
  {
    counters_counter2++;
    timermap.stop();
    return;
  }

  void counters_finalise_()
  {
    std::fstream fstr;
    fstr.open( "counters_log.txt", std::fstream::out | std::fstream::trunc );
    timermap.simpledump( fstr );
    fstr.close();
    return;
  }

}
