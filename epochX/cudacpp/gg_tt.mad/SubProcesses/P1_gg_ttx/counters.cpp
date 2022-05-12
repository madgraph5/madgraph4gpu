#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>

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
  // Now: fortran=-1, cudacpp=0
  // Eventually: fortran=-1, cuda=0, cpp/none=1, cpp/sse4=2, etc...
  constexpr unsigned int nimplC = 2;
  constexpr unsigned int iimplF2C( int iimplF ){ return iimplF + 1; }
  const char* iimplC2TXT( int iimplC )
  {
    const int iimplF = iimplC - 1;
    switch( iimplF )
    {
    case -1: return "Fortran"; break;
    case +0: return "CudaCpp"; break;
    default: assert( false ); break;
    }
  }

#ifdef COUNTERS_USETIMER
  static mgOnGpu::Timer<TIMERTYPE> program_timer;
  static float program_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> matrix1_timer;
  static float matrix1_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> smatrix1_timer;
  static float smatrix1_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> smatrix1multi_timer[nimplC];
  static float smatrix1multi_totaltime[nimplC] = { 0 };
#endif
  static int matrix1_counter = 0;
  static int smatrix1_counter = 0;
  static int smatrix1multi_counter[nimplC] = { 0 };

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

  void counters_smatrix1multi_start_( const int* iimplF, const int* pnevt )
  {
    const unsigned int iimplC = iimplF2C( *iimplF );
    smatrix1multi_counter[iimplC] += *pnevt;
#ifdef COUNTERS_USETIMER
    smatrix1multi_timer[iimplC].Start();
#endif
    return;
  }

  void counters_smatrix1multi_stop_( const int* iimplF )
  {
    const unsigned int iimplC = iimplF2C( *iimplF );
#ifdef COUNTERS_USETIMER
    smatrix1multi_totaltime[iimplC] += smatrix1multi_timer[iimplC].GetDuration();
#endif
    return;
  }

  void counters_finalise_()
  {
    // Write to file
    FILE* f;
    f = fopen( "counters_log.txt", "w" );
#ifdef COUNTERS_USETIMER
    program_totaltime += program_timer.GetDuration();
    fprintf( f, "PROGRAM       : %9.4fs\n", program_totaltime );
    //fprintf( f, "MATRIX1(a)    : %9.4fs for %8d MATRIX1 calls  => throughput is %8.2E calls/s\n", matrix1_totaltime, matrix1_counter, matrix1_counter / matrix1_totaltime );
    //fprintf( f, "MATRIX1(b)    : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", matrix1_totaltime, smatrix1_counter, smatrix1_counter / matrix1_totaltime );
    //fprintf( f, "SMATRIX1      : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", smatrix1_totaltime, smatrix1_counter, smatrix1_counter / smatrix1_totaltime );
    for ( unsigned int iimplC=0; iimplC<nimplC; iimplC++ )
      if ( smatrix1multi_counter[iimplC] > 0 )
        fprintf( f, "SMATRIX1MULTI : %9.4fs for %8d %7s events => throughput is %8.2E events/s\n", smatrix1multi_totaltime[iimplC], smatrix1multi_counter[iimplC], iimplC2TXT( iimplC ), smatrix1multi_counter[iimplC] / smatrix1multi_totaltime[iimplC] );
#else
    //fprintf( f, "MATRIX1       %7s : %8d calls\n", "", matrix1_counter );
    //fprintf( f, "SMATRIX1      %7s : %8d calls\n", "", smatrix1_counter );
    fprintf( f, "SMATRIX1MULTI %7s : %8d events\n", iimplC2TXT( iimplC ), smatrix1multi_counter[iimplC] );
#endif
    fclose( f );
    // Write to stdout
#ifdef COUNTERS_USETIMER
    printf( "PROGRAM       : %9.4fs\n", program_totaltime );
    //printf( "MATRIX1(a)    : %9.4fs for %8d MATRIX1 calls  => throughput is %8.2E calls/s\n", matrix1_totaltime, matrix1_counter, matrix1_counter / matrix1_totaltime );
    //printf( "MATRIX1(b)    : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", matrix1_totaltime, smatrix1_counter, smatrix1_counter / matrix1_totaltime );
    //printf( "SMATRIX1      : %9.4fs for %8d SMATRIX1 calls => throughput is %8.2E calls/s\n", smatrix1_totaltime, smatrix1_counter, smatrix1_counter / smatrix1_totaltime );
    for ( unsigned int iimplC=0; iimplC<nimplC; iimplC++ )
      if ( smatrix1multi_counter[iimplC] > 0 )
        printf( "SMATRIX1MULTI : %9.4fs for %8d %7s events => throughput is %8.2E events/s\n", smatrix1multi_totaltime[iimplC], smatrix1multi_counter[iimplC], iimplC2TXT( iimplC ), smatrix1multi_counter[iimplC] / smatrix1multi_totaltime[iimplC] );
#else
    printf( "MATRIX1       %7s : %8d calls\n", "", matrix1_counter );
    printf( "SMATRIX1      %7s : %8d calls\n", "", smatrix1_counter );
    printf( "SMATRIX1MULTI %7s : %8d events\n", iimplC2TXT( iimplC ), smatrix1multi_counter[iimplC] );
#endif
    return;
  }
}
