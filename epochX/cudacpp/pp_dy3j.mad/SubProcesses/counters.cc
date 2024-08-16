// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#include "timers.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strlen
#include <iostream>
#include <sstream>
#include <string_view>

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
    constexpr int NCOUNTERSMAX = 20;
    static bool disablecounters = false;
    static bool usechronotimers = false;
    // Overall program timer
    static mgOnGpu::ChronoTimer<TIMERTYPE> program_chronotimer;
    static mgOnGpu::RdtscTimer program_rdtsctimer;
    // Individual timers
    static std::string array_tags[NCOUNTERSMAX + 3];
    static mgOnGpu::ChronoTimer<TIMERTYPE> array_chronotimers[NCOUNTERSMAX + 3];
    static mgOnGpu::RdtscTimer array_rdtsctimers[NCOUNTERSMAX + 3];
    static int array_counters[NCOUNTERSMAX + 3] = { 0 };
  }

  void counters_initialise_()
  {
    using namespace counters;
    if( getenv( "CUDACPP_RUNTIME_DISABLECOUNTERS" ) ) disablecounters = true;
    if( getenv( "CUDACPP_RUNTIME_USECHRONOTIMERS" ) ) usechronotimers = true;
    for( int icounter = 1; icounter < NCOUNTERSMAX + 1; icounter++ )
      array_tags[icounter] = ""; // ensure that this is initialized to ""
    if( usechronotimers ) program_chronotimer.start();
    else program_rdtsctimer.start();
    return;
  }

  void counters_register_counter_( const int* picounter, const char* ctag )
  {
    using namespace counters;
    int icounter = *picounter;
    std::cout << "INFO: register counter #" << icounter << " with tag '" << ctag << "' (tag strlen=" << strlen( ctag ) << ")" << std::endl;
    const std::string tag( ctag );
    if( icounter < 1 || icounter >= NCOUNTERSMAX + 1 )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Invalid counter# '" << icounter << "' (valid values are 1 to " << NCOUNTERSMAX << ")";
      throw std::runtime_error( sstr.str() );
    }
    if( tag == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Invalid empty tag ''";
      throw std::runtime_error( sstr.str() );
    }
    if( array_tags[icounter] == "" )
    {
      array_tags[icounter] = tag;
    }
    else
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " already exists with tag '" << array_tags[icounter] << "'";
      throw std::runtime_error( sstr.str() );
    }
    return;
  }

  void counters_start_counter_( const int* picounter, const int* pnevt )
  {
    using namespace counters;
    if( disablecounters ) return;
    int icounter = *picounter;
    if( array_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    array_counters[icounter] += *pnevt;
    if( usechronotimers ) array_chronotimers[icounter].start();
    else array_rdtsctimers[icounter].start();
    return;
  }

  void counters_stop_counter_( const int* picounter )
  {
    using namespace counters;
    if( disablecounters ) return;
    int icounter = *picounter;
    if( array_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    if( usechronotimers ) array_chronotimers[icounter].stop();
    else array_rdtsctimers[icounter].stop();
    return;
  }

  inline bool starts_with( std::string_view str, std::string_view prefix ) // https://stackoverflow.com/a/42844629
  {
    return str.size() >= prefix.size() && str.compare( 0, prefix.size(), prefix ) == 0;
  }

  inline bool ends_with( std::string_view str, std::string_view suffix ) // https://stackoverflow.com/a/42844629
  {
    return str.size() >= suffix.size() && str.compare( str.size() - suffix.size(), suffix.size(), suffix ) == 0;
  }

  void counters_finalise_()
  {
    using namespace counters;
    // Dump program counters
    if( usechronotimers ) program_chronotimer.stop();
    else program_rdtsctimer.stop();
    float program_totaltime = ( usechronotimers ? program_chronotimer.getDurationSeconds() : program_rdtsctimer.getDurationSeconds() );
    if( usechronotimers ) printf( " [COUNTERS] *** USING STD::CHRONO TIMERS ***\n" );
    else printf( " [COUNTERS] *** USING RDTSC-BASED TIMERS ***\n" );
    printf( " [COUNTERS] PROGRAM TOTAL                         : %9.4fs\n", program_totaltime );
    if( disablecounters ) return;
    // Extract time duration from all timers
    float array_totaltimes[NCOUNTERSMAX + 3] = { 0 };
    for( int icounter = 1; icounter < NCOUNTERSMAX + 1; icounter++ )
    {
      if( usechronotimers )
        array_totaltimes[icounter] = array_chronotimers[icounter].getDurationSeconds();
      else
        array_totaltimes[icounter] = array_rdtsctimers[icounter].getDurationSeconds();
    }
    // Create counter[0] "Fortran Other"
    array_tags[0] = "Fortran Other";
    array_counters[0] = 1;
    array_totaltimes[0] = program_totaltime;
    for( int icounter = 1; icounter < NCOUNTERSMAX + 1; icounter++ )
    {
      if( !starts_with( array_tags[icounter], "PROGRAM" ) ) // skip counters whose tags start with "PROGRAM"
        array_totaltimes[0] -= array_totaltimes[icounter];
    }
    // Create counters[NCOUNTERSMAX+2] "OVERALL MEs" and counters[NCOUNTERSMAX+1] "OVERALL NON-MEs"
    array_tags[NCOUNTERSMAX + 2] = "OVERALL MEs";
    array_counters[NCOUNTERSMAX + 2] = 0;
    array_totaltimes[NCOUNTERSMAX + 2] = 0;
    for( int icounter = 1; icounter < NCOUNTERSMAX + 1; icounter++ )
    {
      if( ends_with( array_tags[icounter], "MEs" ) ) // include counters whose tags end with "MEs"
      {
        array_counters[NCOUNTERSMAX + 2] += array_counters[icounter];
        array_totaltimes[NCOUNTERSMAX + 2] += array_totaltimes[icounter];
      }
    }
    array_tags[NCOUNTERSMAX + 1] = "OVERALL NON-MEs";
    array_counters[NCOUNTERSMAX + 1] = 1;
    array_totaltimes[NCOUNTERSMAX + 1] = program_totaltime - array_totaltimes[NCOUNTERSMAX + 2];
    // Dump individual counters
    for( int icounter = 0; icounter < NCOUNTERSMAX + 3; icounter++ )
    {
      if( array_tags[icounter] != "" )
      {
        if( array_counters[icounter] > 1 ) // event counters
        {
          printf( " [COUNTERS] %-30s ( %2d ) : %9.4fs for %8d events => throughput is %8.2E events/s\n",
                  array_tags[icounter].c_str(),
                  icounter,
                  array_totaltimes[icounter],
                  array_counters[icounter],
                  array_totaltimes[icounter] / array_counters[icounter] );
        }
        else if( array_counters[icounter] == 1 ) // one-off counters for initialisation tasks (e.g. helicity filtering)
        {
          printf( " [COUNTERS] %-30s ( %2d ) : %9.4fs\n",
                  array_tags[icounter].c_str(),
                  icounter,
                  array_totaltimes[icounter] );
        }
      }
    }
    return;
  }
}
