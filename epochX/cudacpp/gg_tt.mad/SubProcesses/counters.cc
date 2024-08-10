// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>
#include <cstring> // for strlen
#include <map>
#include <memory>

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
    // Overall program timer
    static mgOnGpu::Timer<TIMERTYPE> program_timer;
    static float program_totaltime = 0;
    // Individual timers
    static std::map<unsigned int, std::string> map_tags;
    static std::map<unsigned int, std::unique_ptr< mgOnGpu::Timer<TIMERTYPE> > > map_timers;
    static std::map<unsigned int, float > map_totaltimes;
    static std::map<unsigned int, int > map_counters;
  }
  
  void counters_initialise_()
  {
    using namespace counters;
    program_timer.Start();
    return;
  }

  void counters_register_counter_( const int* picounter, const char* ctag )
  {
    using namespace counters;
    unsigned int icounter = *picounter;
    std::cout << "INFO: register counter #" << icounter << " with tag '" << ctag << "' (tag strlen=" << strlen(ctag) << ")" << std::endl;
    const std::string tag(ctag);
    if( map_tags.find( icounter ) == map_tags.end() )
    {
      map_tags[icounter] = tag;
      map_timers[icounter] = std::make_unique<mgOnGpu::Timer<TIMERTYPE>>();
      map_totaltimes[icounter] = 0;
      map_counters[icounter] = 0;
    }
    else
    {
      std::cout << "ERROR! counter #" << icounter << " already exists with tag '" << map_tags[ icounter ] << "'" << std::endl;
    }
    return;
  }

  void counters_start_counter_( const int* picounter, const int* pnevt )
  {
    using namespace counters;
    unsigned int icounter = *picounter;
    if( map_tags.find( icounter ) == map_tags.end() )
      std::cout << "ERROR! counter #" << icounter << " does not exist" << std::endl;
    map_counters[icounter] += *pnevt;
    map_timers[icounter]->Start();
    return;
  }

  void counters_stop_counter_( const int* picounter )
  {
    using namespace counters;
    unsigned int icounter = *picounter;
    if( map_tags.find( icounter ) == map_tags.end() )
      std::cout << "ERROR! counter #" << icounter << " does not exist" << std::endl;
    map_totaltimes[icounter] += map_timers[icounter]->GetDuration();
    return;
  }

  void counters_finalise_()
  {
    using namespace counters;
    program_totaltime += program_timer.GetDuration();
    // Write to stdout
    float overhead_totaltime = program_totaltime;
    for( auto const& [icounter, totaltime] : map_totaltimes ) overhead_totaltime -= totaltime;
    printf( " [COUNTERS] PROGRAM TOTAL          : %9.4fs\n", program_totaltime );
    printf( " [COUNTERS] Fortran Overhead ( 0 ) : %9.4fs\n", overhead_totaltime );
    for( auto const& [icounter, tag] : map_tags )
    {
      if( map_counters[icounter] > 1 ) // event counters
      {
        printf( " [COUNTERS] %11s      ( %1d ) : %9.4fs for %8d events => throughput is %8.2E events/s\n",
                tag.c_str(),
                icounter,
                map_totaltimes[icounter],
                map_counters[icounter],
                map_totaltimes[icounter] / map_counters[icounter] );
      }
      else if( map_counters[icounter] == 1 ) // one-off counters for initialisation tasks (e.g. helicity filtering)
      {
        printf( " [COUNTERS] %11s      ( %1d ) : %9.4fs\n",
                tag.c_str(),
                icounter,
                map_totaltimes[icounter] );
      }
    }
    return;
  }
}
