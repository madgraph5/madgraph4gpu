// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>
#include <cstring> // for strlen
#include <sstream>

// NB1: The C functions counters_xxx_ in this file are called by Fortran code
// Hence the trailing "_": 'call counters_end()' links to counters_end_
// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

// NB2: This file also contains C++ code and is built using g++
// Hence use 'extern "C"' to avoid name mangling by the C++ compiler
// See https://www.geeksforgeeks.org/extern-c-in-c

extern "C"
{
  namespace counters
  {
    constexpr int NCOUNTERSMAX = 10;
    // Overall program timer
    static mgOnGpu::Timer<TIMERTYPE> program_timer;
    static float program_totaltime = 0;
    // Individual timers
    static std::string map_tags[NCOUNTERSMAX];
    static mgOnGpu::Timer<TIMERTYPE> map_timers[NCOUNTERSMAX];
    static float map_totaltimes[NCOUNTERSMAX] = { 0 };
    static int map_counters[NCOUNTERSMAX] = { 0 };
  }
  
  void counters_initialise_()
  {
    using namespace counters;
    for( int icounter=0; icounter<NCOUNTERSMAX; icounter++ )
      map_tags[icounter] = ""; // ensure that this is initialized to ""
    program_timer.Start();
    return;
  }

  void counters_register_counter_( const int* picounter, const char* ctag )
  {
    using namespace counters;
    int icounter = *picounter;
    std::cout << "INFO: register counter #" << icounter << " with tag '" << ctag << "' (tag strlen=" << strlen(ctag) << ")" << std::endl;
    const std::string tag(ctag);
    if( icounter < 0 || icounter >= NCOUNTERSMAX )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Invalid counter# '" << icounter << "'";
      throw std::runtime_error( sstr.str() );
    }
    if( tag == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Invalid empty tag ''";
      throw std::runtime_error( sstr.str() );
    }
    if( map_tags[icounter] == "" )
    {
      map_tags[icounter] = tag;
    }
    else
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " already exists with tag '" << map_tags[ icounter ] << "'";
      throw std::runtime_error( sstr.str() );
    }
    return;
  }

  void counters_start_counter_( const int* picounter, const int* pnevt )
  {
    using namespace counters;
    int icounter = *picounter;
    if( map_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    map_counters[icounter] += *pnevt;
    map_timers[icounter].Start();
    return;
  }

  void counters_stop_counter_( const int* picounter )
  {
    using namespace counters;
    int icounter = *picounter;
    if( map_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    map_totaltimes[icounter] += map_timers[icounter].GetDuration();
    return;
  }

  void counters_finalise_()
  {
    using namespace counters;
    program_totaltime += program_timer.GetDuration();
    // Write to stdout
    float overhead_totaltime = program_totaltime;
    for( int icounter=0; icounter<NCOUNTERSMAX; icounter++ ) overhead_totaltime -= map_totaltimes[icounter];
    printf( " [COUNTERS] PROGRAM TOTAL          : %9.4fs\n", program_totaltime );
    printf( " [COUNTERS] Fortran Overhead ( 0 ) : %9.4fs\n", overhead_totaltime );
    for( int icounter=0; icounter<NCOUNTERSMAX; icounter++ )
    {
      if( map_tags[icounter] != "" )
      {
        if( map_counters[icounter] > 1 ) // event counters
        {
          printf( " [COUNTERS] %11s      ( %1d ) : %9.4fs for %8d events => throughput is %8.2E events/s\n",
                  map_tags[icounter].c_str(),
                  icounter,
                  map_totaltimes[icounter],
                  map_counters[icounter],
                  map_totaltimes[icounter] / map_counters[icounter] );
        }
        else if( map_counters[icounter] == 1 ) // one-off counters for initialisation tasks (e.g. helicity filtering)
        {
          printf( " [COUNTERS] %11s      ( %1d ) : %9.4fs\n",
                  map_tags[icounter].c_str(),
                  icounter,
                  map_totaltimes[icounter] );
        }
      }
    }
    return;
  }
}
