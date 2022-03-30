#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <stdio.h>

#undef COUNTERS_USETIMER
#define COUNTERS_USETIMER 1 // comment out to disable timers and keep only counters

// NB1: The C functions counters_xxx_ in this file are called by Fortran code
// Hence the trailing "_": 'call counters_end()' links to counters_end_
// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

// NB2: This file also contains C++ code and is built using g++
// Hence use 'extern "C"' to avoid name mangling by the C++ compiler
// See https://www.geeksforgeeks.org/extern-c-in-c

extern "C"
{
#ifdef COUNTERS_USETIMER
  static mgOnGpu::Timer<TIMERTYPE> program_timer;
  static float program_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> matrix1_timer;
  static float matrix1_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> smatrix1_timer;
  static float smatrix1_totaltime = 0;
#endif
  static int matrix1_counter = 0;
  static int smatrix1_counter = 0;

  void counters_initialise_()
  {
#ifdef COUNTERS_USETIMER
    program_timer.Start();
#endif
    return;
  }

  void counters_matrix1_start_()
  {
    matrix1_counter++;
#ifdef COUNTERS_USETIMER
    matrix1_timer.Start();
#endif
    return;
  }

  void counters_matrix1_stop_()
  {
#ifdef COUNTERS_USETIMER
    matrix1_totaltime += matrix1_timer.GetDuration();
#endif
    return;
  }

  void counters_smatrix1_start_()
  {
    smatrix1_counter++;
#ifdef COUNTERS_USETIMER
    smatrix1_timer.Start();
#endif
    return;
  }

  void counters_smatrix1_stop_()
  {
#ifdef COUNTERS_USETIMER
    smatrix1_totaltime += smatrix1_timer.GetDuration();
#endif
    return;
  }

  void counters_finalise_()
  {
    FILE* f;
    f = fopen( "counters_log.txt", "w" );
#ifdef COUNTERS_USETIMER
    program_totaltime += program_timer.GetDuration();
    fprintf( f, "PROGRAM    : %9.4fs\n", program_totaltime );
    fprintf( f, "MATRIX1(a) : %9.4fs for %8d MATRIX1 calls  => throughput is %8.2E calls/s\n", matrix1_totaltime, matrix1_counter, matrix1_counter / matrix1_totaltime );
    fprintf( f, "MATRIX1(b) : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", matrix1_totaltime, smatrix1_counter, smatrix1_counter / matrix1_totaltime );
    fprintf( f, "SMATRIX1   : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", smatrix1_totaltime, smatrix1_counter, smatrix1_counter / smatrix1_totaltime );
#else
    fprintf( f, "MATRIX1  : %8d calls\n", matrix1_counter );
    fprintf( f, "SMATRIX1 : %8d calls\n", matrix1_counter );
#endif
    fclose( f );
    return;
  }
}
