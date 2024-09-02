// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: S. Hageboeck, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strlen
#include <iostream>
#include <set>
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
    constexpr int NCOUNTERSMAX = 30;
    static bool disablecalltimers = false;
    static bool disabletesttimers = false;
#ifdef MGONGPU_HASRDTSC
    static bool usechronotimers = false;
#else
    constexpr bool usechronotimers = true;
#endif
    static bool removetimeroverhead = false;
    // Overall program timer
    static mgOnGpu::ChronoTimer<TIMERTYPE> program_chronotimer;
    static mgOnGpu::RdtscTimer program_rdtsctimer;
    // Individual timers
    static std::string array_tags[NCOUNTERSMAX + 4];
    static bool array_istesttimer[NCOUNTERSMAX + 4];
    static mgOnGpu::ChronoTimer<TIMERTYPE> array_chronotimers[NCOUNTERSMAX + 4];
    static mgOnGpu::RdtscTimer array_rdtsctimers[NCOUNTERSMAX + 4];
    static int array_counters[NCOUNTERSMAX + 4] = { 0 };
    static std::set<int> array_included[NCOUNTERSMAX + 4];
    static int array_stopcalls[NCOUNTERSMAX + 4] = { 0 };
    static float overheadpercallseconds = 0;
  }

  inline bool starts_with( std::string_view str, std::string_view prefix ) // https://stackoverflow.com/a/42844629
  {
    return str.size() >= prefix.size() && str.compare( 0, prefix.size(), prefix ) == 0;
  }

  inline bool ends_with( std::string_view str, std::string_view suffix ) // https://stackoverflow.com/a/42844629
  {
    return str.size() >= suffix.size() && str.compare( str.size() - suffix.size(), suffix.size(), suffix ) == 0;
  }

  void counters_register_counter_( const int* picounter, const char* ctag )
  {
    using namespace counters;
    int icounter = *picounter;
    std::cout << "INFO: register counter #" << icounter << " with tag '" << ctag << "' (tag strlen=" << strlen( ctag ) << ")" << std::endl;
    const std::string tag( ctag );
    if ( icounter < 1 || icounter == NCOUNTERSMAX + 1 || icounter == NCOUNTERSMAX + 2 || icounter >= NCOUNTERSMAX + 4 ) // allow icalibcounter = NCOUNTERSMAX+3
    {
      std::ostringstream sstr;
      sstr << "ERROR! Invalid counter # '" << icounter << "' (valid values are 1 to " << NCOUNTERSMAX << ")";
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
      if( starts_with( array_tags[icounter], "TEST" ) ) array_istesttimer[icounter] = true;
    }
    else
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " already exists with tag '" << array_tags[icounter] << "'";
      throw std::runtime_error( sstr.str() );
    }
    return;
  }

  void counters_include_counter1_within_counter2_( const int* picounter1, const int* picounter2 )
  {
    using namespace counters;
    int icounter1 = *picounter1;
    int icounter2 = *picounter2;
    std::cout << "INFO: declare that counter #" << icounter1 << " is included within counter #" << icounter2 << std::endl;
    if( array_tags[icounter1] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Counter # '" << icounter1 << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    if( array_tags[icounter2] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Counter # '" << icounter2 << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    // For the moment only the inclusion of a TEST timer into a non-TEST timer is allowed
    // (the only aim of this mechanism is the removal of overheads, which otherwise would be too complex)
    if( !array_istesttimer[icounter1] )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Counter # '" << icounter1 << " is not a TEST timer";
      throw std::runtime_error( sstr.str() );
    }
    if( array_istesttimer[icounter2] )
    {
      std::ostringstream sstr;
      sstr << "ERROR! Counter # '" << icounter2 << " is a TEST timer";
      throw std::runtime_error( sstr.str() );
    }
    array_included[icounter2].insert( icounter1 );
    return;
  }

  void counters_start_counter_( const int* picounter, const int* pnevt )
  {
    using namespace counters;
    if( disablecalltimers ) return;
    int icounter = *picounter;
    if( disabletesttimers && array_istesttimer[icounter] ) return;
    if( array_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    array_counters[icounter] += *pnevt;
    if( usechronotimers )
      array_chronotimers[icounter].start();
    else
      array_rdtsctimers[icounter].start();
    return;
  }

  void counters_stop_counter_( const int* picounter )
  {
    using namespace counters;
    if( disablecalltimers ) return;
    int icounter = *picounter;
    if( disabletesttimers && array_istesttimer[icounter] ) return;
    if( array_tags[icounter] == "" )
    {
      std::ostringstream sstr;
      sstr << "ERROR! counter #" << icounter << " does not exist";
      throw std::runtime_error( sstr.str() );
    }
    if( usechronotimers )
      array_chronotimers[icounter].stop();
    else
      array_rdtsctimers[icounter].stop();
    array_stopcalls[icounter]++;
    return;
  }

  void counters_initialise_()
  {
    using namespace counters;
#ifdef MGONGPU_HASRDTSC
    if( getenv( "CUDACPP_RUNTIME_USECHRONOTIMERS" ) ) usechronotimers = true;
#endif
    if( getenv( "CUDACPP_RUNTIME_REMOVECOUNTEROVERHEAD" ) ) removetimeroverhead = true;
    for( int icounter = 0; icounter < NCOUNTERSMAX + 4; icounter++ ) // include icalibcounter = NCOUNTERSMAX+3
    {
      array_tags[icounter] = "";           // ensure that this is initialized to ""
      array_istesttimer[icounter] = false; // ensure that this is initialized to false
    }
    if( usechronotimers )
      program_chronotimer.start();
    else
      program_rdtsctimer.start();
    if( removetimeroverhead )
    {
      const int icalibcounter = NCOUNTERSMAX + 3;
      const int nevtdummy = 1;
      counters_register_counter_( &icalibcounter, "OVERHEAD CALIBRATION" );
      mgOnGpu::ChronoTimer<std::chrono::high_resolution_clock> calibtimer;
      calibtimer.start();
      constexpr size_t ncall = 10000000; // 10M calls are expected to take slightly less than ~1s (this will be in counter overhead)
      for( size_t icall = 0; icall < ncall; icall++ )
      {
        counters_start_counter_( &icalibcounter, &nevtdummy );
        counters_stop_counter_( &icalibcounter );
      }
      calibtimer.stop();
      overheadpercallseconds = calibtimer.getTotalDurationSeconds() / ncall;
    }
    if( getenv( "CUDACPP_RUNTIME_DISABLECALLTIMERS" ) ) disablecalltimers = true;
    if( getenv( "CUDACPP_RUNTIME_DISABLETESTTIMERS" ) ) disabletesttimers = true;
    return;
  }

  void counters_finalise_()
  {
    using namespace counters;
    // Compute program counters
    if( usechronotimers )
      program_chronotimer.stop();
    else
      program_rdtsctimer.stop();
    float program_totaltime = ( usechronotimers ? program_chronotimer.getTotalDurationSeconds() : program_rdtsctimer.getTotalDurationSeconds() );
    float program_overhead = 0;
    // Extract time duration from all timers
    float array_totaltimes[NCOUNTERSMAX + 4] = { 0 };
    float array_overheads[NCOUNTERSMAX + 4] = { 0 };
    for( int icounter = 1; icounter < NCOUNTERSMAX + 4; icounter++ ) // include icalibcounter = NCOUNTERSMAX+3
    {
      if( icounter == NCOUNTERSMAX + 1 ) continue;
      if( icounter == NCOUNTERSMAX + 2 ) continue;
      if( usechronotimers )
        array_totaltimes[icounter] = array_chronotimers[icounter].getTotalDurationSeconds();
      else
        array_totaltimes[icounter] = array_rdtsctimers[icounter].getTotalDurationSeconds();
      // Remove own overheads if required
      if( removetimeroverhead )
      {
        array_overheads[icounter] = array_stopcalls[icounter] * overheadpercallseconds;
        array_totaltimes[icounter] -= array_overheads[icounter];
        program_overhead += array_overheads[icounter]; // remove overheads of all timers including TEST timers
      }
    }
    // Remove overheads of included timers if any
    if( removetimeroverhead )
    {
      for( int icounter = 1; icounter < NCOUNTERSMAX + 1; icounter++ ) // no need to include icalibcounter = NCOUNTERSMAX+3
      {
        for( int icounterIn : array_included[icounter] )
          array_totaltimes[icounter] -= array_overheads[icounterIn];
      }
      printf( " INFO: COUNTERS overhead : %9.4fs for 1M start/stop cycles\n", overheadpercallseconds * 1E6 );
    }
    // Dump program counters (after removing overhead if required)
    program_totaltime -= program_overhead;
    std::string timertypemsg = ( usechronotimers ? "STD::CHRONO" : "RDTSC-BASED" );
    std::string overheadmsg = ( removetimeroverhead ? "(remove timer overhead)" : "(do not remove timer overhead)" );
    if( removetimeroverhead )
    {
      printf( " [COUNTERS] PROGRAM TOTAL+COUNTEROVERHEAD         : %9.4fs\n", program_totaltime + program_overhead );
      printf( " [COUNTERS] PROGRAM COUNTEROVERHEAD               : %9.4fs\n", program_overhead );
      printf( " -------------------------------------------------------------\n" );
    }    
    printf( " [COUNTERS] *** USING %s TIMERS %s ***\n", timertypemsg.c_str(), overheadmsg.c_str() );
    printf( " [COUNTERS] PROGRAM TOTAL                         : %9.4fs\n", program_totaltime );
    if( disablecalltimers ) return;
    // Create counter[0] "Fortran Other"
    array_tags[0] = "Fortran Other";
    array_counters[0] = 1;
    array_totaltimes[0] = program_totaltime;
    for( int icounter = 1; icounter < NCOUNTERSMAX + 4; icounter++ ) // include icalibcounter = NCOUNTERSMAX+3
    {
      if( !array_istesttimer[icounter] ) // skip TEST counters
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
    for( int icounter = 0; icounter < NCOUNTERSMAX + 3; icounter++ ) // exclude icalibcounter = NCOUNTERSMAX+3 (would print a negative value here!)
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
                  array_counters[icounter] / array_totaltimes[icounter] );
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
