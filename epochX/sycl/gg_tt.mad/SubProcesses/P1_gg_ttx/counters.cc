// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.
//
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2022-2023) for the MG5aMC SYCL plugin.

#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>

// NB1: The C functions counters_xxx_ in this file are called by Fortran code
// Hence the trailing "_": 'call counters_end()' links to counters_end_
// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

// NB2: This file also contains C++ code and is built using g++
// Hence use 'extern "C"' to avoid name mangling by the C++ compiler
// See https://www.geeksforgeeks.org/extern-c-in-c

extern "C"
{
  // Now: fortran=-1, sycl=0
  // Eventually: fortran=-1, sycl=0, cpp/none=1, cpp/sse4=2, etc...
  constexpr unsigned int nimplC = 2;
  constexpr unsigned int iimplF2C( int iimplF ){ return iimplF + 1; }
  const char* iimplC2TXT( int iimplC )
  {
    const int iimplF = iimplC - 1;
    switch( iimplF )
    {
    case -1: return "Fortran"; break;
    case +0: return "SYCL"; break;
    default: assert( false ); break;
    }
  }

  static mgOnGpu::Timer<TIMERTYPE> program_timer;
  static float program_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> matrix1_timer;
  static float matrix1_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> smatrix1_timer;
  static float smatrix1_totaltime = 0;
  static mgOnGpu::Timer<TIMERTYPE> smatrix1multi_timer[nimplC];
  static float smatrix1multi_totaltime[nimplC] = { 0 };
  static int matrix1_counter = 0;
  static int smatrix1_counter = 0;
  static int smatrix1multi_counter[nimplC] = { 0 };

  void counters_initialise_()
  {
    program_timer.Start();
    return;
  }

  void counters_matrix1_start_()
  {
    matrix1_counter++;
    matrix1_timer.Start();
    return;
  }

  void counters_matrix1_stop_()
  {
    matrix1_totaltime += matrix1_timer.GetDuration();
    return;
  }

  void counters_smatrix1_start_()
  {
    smatrix1_counter++;
    smatrix1_timer.Start();
    return;
  }

  void counters_smatrix1_stop_()
  {
    smatrix1_totaltime += smatrix1_timer.GetDuration();
    return;
  }

  void counters_smatrix1multi_start_( const int* iimplF, const int* pnevt )
  {
    const unsigned int iimplC = iimplF2C( *iimplF );
    smatrix1multi_counter[iimplC] += *pnevt;
    smatrix1multi_timer[iimplC].Start();
    return;
  }

  void counters_smatrix1multi_stop_( const int* iimplF )
  {
    const unsigned int iimplC = iimplF2C( *iimplF );
    smatrix1multi_totaltime[iimplC] += smatrix1multi_timer[iimplC].GetDuration();
    return;
  }

  void counters_finalise_()
  {
    program_totaltime += program_timer.GetDuration();
    // Write to stdout
    float overhead_totaltime = program_totaltime;
    for ( unsigned int iimplC=0; iimplC<nimplC; iimplC++ ) overhead_totaltime -= smatrix1multi_totaltime[iimplC];
    printf( " [COUNTERS] PROGRAM TOTAL          : %9.4fs\n", program_totaltime );
    printf( " [COUNTERS] Fortran Overhead ( 0 ) : %9.4fs\n", overhead_totaltime );
    for ( unsigned int iimplC=0; iimplC<nimplC; iimplC++ )
      if ( smatrix1multi_counter[iimplC] > 0 )
        printf( " [COUNTERS] %7s MEs      ( %1d ) : %9.4fs for %8d events => throughput is %8.2E events/s\n",
                iimplC2TXT( iimplC ),
                iimplC + 1,
                smatrix1multi_totaltime[iimplC],
                smatrix1multi_counter[iimplC],
                smatrix1multi_counter[iimplC] / smatrix1multi_totaltime[iimplC] );
    return;
  }
}
