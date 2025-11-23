// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jul 2020) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2025) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUTIMERMAP2_H
#define MGONGPUTIMERMAP2_H 1

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
//#include "nvtx.h"
//#pragma GCC diagnostic pop

#include "timer2.h"
#define TIMERTYPE std::chrono::high_resolution_clock

namespace mgOnGpu
{
  class TimerMap2
  {

  public:

    // Constructor
    TimerMap2()
      : m_chronoTimer()
      , m_rdtscTimer()
      , m_partitionIdToKey()
      , m_active( 0 )
      , m_partitionTotalCounts()
      , m_useChronoTimers( false )
      , m_started( false )
    {
#ifdef MGONGPU_HASRDTSC
      if( getenv( "CUDACPP_RUNTIME_USECHRONOTIMERS" ) ) m_useChronoTimers = true;
#else
      m_useChronoTimers = true;
#endif
    }

    // Destructor
    virtual ~TimerMap2() {}

    // Add a partition
    void addPartition( size_t id, const std::string& key )
    {
      assert( id > 0 ); // id == 0 signals that no partition is active
      assert( m_partitionIdToKey.find( id ) == m_partitionIdToKey.end() );
      for( auto ip: m_partitionIdToKey ) assert( ip.second != key );
      m_partitionIdToKey[id] = key;
      m_partitionTotalCounts[id] = 0;
    }

    // Start the timer for a specific partition (key must be a non-empty string)
    // Stop the timer for the current partition if there is one active
    uint64_t start( size_t id )
    {
      assert( id > 0 );
      //assert( m_partitionIdToKey.find( id ) != m_partitionIdToKey.end() ); // unnecessary overhead
      // Close the previously active partition
      uint64_t last = stop();
      // Switch to a new partition
      if( !m_started )
      {
        if( m_useChronoTimers )
          m_chronoTimer.start();
        else
          m_rdtscTimer.start();
        m_started = true;
      }
      m_active = id;
      // Open a new Cuda NVTX range
      //NVTX_PUSH( m_partitionIdToKey[id].c_str(), id ); // unnecessary overhead
      // Return last duration
      return last;
    }

    // Stop the timer for the current partition if there is one active
    uint64_t stop()
    {
      // Close the previously active partition
      uint64_t last = 0;
      if( m_started )
      {
        if( m_useChronoTimers )
          last = m_chronoTimer.getCountsSinceStart();
        else
          last = m_rdtscTimer.getCountsSinceStart();
        m_partitionTotalCounts[m_active] += last;
        if( m_useChronoTimers )
          m_chronoTimer.stop();
        else
          m_rdtscTimer.stop();
        m_started = false;
      }
      m_active = 0;
      // Close the current Cuda NVTX range
      //NVTX_POP(); // unnecessary overhead
      // Return last duration
      return last;
    }

    // Return timer calibration (at this point in time for rdtsc, constant in time for chrono)
    float secondsPerCount()
    {
      if( m_useChronoTimers )
        return m_chronoTimer.secondsPerCount();
      else
        return m_rdtscTimer.secondsPerCount();
    }

    // Dump the overall results
    void dump( const std::string totalKey = "TOTAL", std::ostream& ostr = std::cout )
    {
      // Improve key formatting
      size_t maxsize = 0;
      for( auto ip: m_partitionIdToKey )
        maxsize = std::max( maxsize, ip.second.size() );
      maxsize = std::max( maxsize, totalKey.size() );
      // Compute individual partition total times from partition total counts
      std::map<std::string, float> partitionTotalTimes;
      float secPerCount = secondsPerCount();
      for( auto ip: m_partitionTotalCounts )
      {
        std::string key = m_partitionIdToKey[ip.first];
        partitionTotalTimes[key] = m_partitionTotalCounts[ip.first] * secPerCount;
      }
      // Compute the overall total
      float total = 0;
      for( auto ip: partitionTotalTimes ) total += ip.second;
      // Dump individual partition timers and the overall total
      // NB: 'setw' affects only the next field (of any type)
      ostr << std::setprecision( 6 ); // set precision (default=6): affects all floats
      ostr << std::fixed;             // fixed format: affects all floats
      for( auto ip: partitionTotalTimes )
        ostr << std::setw( maxsize ) << ip.first << " : "
             << std::setw( 12 ) << ip.second << " sec" << std::endl;
      ostr << std::setw( maxsize ) << totalKey << " : "
           << std::setw( 12 ) << total << " sec" << std::endl;
      ostr << std::defaultfloat; // default format: affects all floats
    }

  private:

    ChronoTimer<TIMERTYPE> m_chronoTimer;
    RdtscTimer m_rdtscTimer;
    std::map<size_t, std::string> m_partitionIdToKey;
    size_t m_active;
    std::map<size_t, uint64_t> m_partitionTotalCounts;
    bool m_useChronoTimers;
    bool m_started; // when the timer is stopped, it must be explicitly restarted
  };

}

#endif // MGONGPUTIMERMAP2_H
