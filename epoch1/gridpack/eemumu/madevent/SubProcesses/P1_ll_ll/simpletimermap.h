#ifndef MGONGPUSIMPLETIMERMAP_H
#define MGONGPUSIMPLETIMERMAP_H 1

#include <cassert>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>

#ifdef USE_NVTX
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include "nvtx.h"
#pragma GCC diagnostic pop
#endif

#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

namespace mgOnGpu
{
  class SimpleTimerMap
  {

  public:

    SimpleTimerMap() : m_timer(), m_active(""), m_partitionIds(), m_partitionCounters(), m_partitionTimers() {}
    virtual ~SimpleTimerMap() {}

    // Start the timer for a specific partition (key must be a non-empty string)
    // Stop the timer for the current partition if there is one active
    float start( const std::string& key )
    {
      assert( key != "" );
      // Close the previously active partition
      float last = stop();
      // Switch to a new partition
      m_timer.Start();
      m_active = key;
      if( m_partitionIds.find(key) == m_partitionIds.end() )
      {
        m_partitionIds[key] = m_partitionTimers.size();
        m_partitionCounters[key] = 0;
        m_partitionTimers[key] = 0;
      }
      m_partitionCounters[key]++;
#ifdef USE_NVTX
      // Open a new Cuda NVTX range
      NVTX_PUSH( key.c_str(), m_partitionIds[key] );
#endif
      // Return last duration
      return last;
    }

    // Stop the timer for the current partition if there is one active
    float stop()
    {
      // Close the previously active partition
      float last = 0;
      if ( m_active != "" )
      {
        last = m_timer.GetDuration();
        m_partitionTimers[m_active] += last;
      }
      m_active = "";
#ifdef USE_NVTX
      // Close the current Cuda NVTX range
      NVTX_POP();
#endif
      // Return last duration
      return last;
    }

    // Dump the overall results
    void simpledump( std::ostream& ostr = std::cout )
    {
      // Improve key formatting
      size_t maxsize = 0;
      for ( auto ip : m_partitionTimers )
        maxsize = std::max( maxsize, ip.first.size() );
      // Dump individual partition timers
      // NB: 'setw' affects only the next field (of any type)
      ostr << std::setprecision(4); // set precision (default=6): affects all floats
      ostr << std::fixed; // fixed format: affects all floats
      for ( auto ip : m_partitionTimers )
        ostr << std::setw(maxsize) << ip.first << " : "
             << std::setw(8) << m_partitionCounters[ip.first] << " calls in "
             << std::setw(8) << ip.second << " sec"
             << " : throughput " << std::scientific
             << m_partitionCounters[ip.first] / ip.second << " calls/sec" << std::endl;
#if __GNUC__ > 4
      ostr << std::defaultfloat; // default format: affects all floats
#endif
    }

    private:

    Timer<TIMERTYPE> m_timer;
    std::string m_active;
    std::map< std::string, uint32_t > m_partitionIds;
    std::map< std::string, uint32_t > m_partitionCounters;
    std::map< std::string, float > m_partitionTimers;

    };

    }

#endif // MGONGPUSIMPLETIMERMAP_H
