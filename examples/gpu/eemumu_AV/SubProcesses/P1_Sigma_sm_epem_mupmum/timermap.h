#ifndef MGONGPUTIMERMAP_H 
#define MGONGPUTIMERMAP_H 1

#include <cassert>
#include <map>
#include <string>

#include "nvtx.h"
#include "timer.h"
#define TIMERTYPE std::chrono::high_resolution_clock

namespace mgOnGpu
{
  class TimerMap
  {

  public:

    TimerMap() : m_timer(), m_active(""), m_partitionTimers(), m_partitionIds() {}
    virtual ~TimerMap() {}

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
      if( m_partitionTimers.find(key) == m_partitionTimers.end() )
      {
        m_partitionIds[key] = m_partitionTimers.size();
        m_partitionTimers[key] = 0;
      }
      // Open a new Cuda NVTX range
      NVTX_PUSH( key.c_str(), m_partitionIds[key] );
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
      // Close the current Cuda NVTX range
      NVTX_POP();
      // Return last duration
      return last;
    }

    // Dump the overall results
    void dump()
    {
      // Improve key formatting
      const std::string totalKey = "TOTAL      ";
      const std::string total2Key = "TOTAL(n-2) ";
      size_t maxsize = 0;
      for ( auto ip : m_partitionTimers )
        maxsize = std::max( maxsize, ip.first.size() );
      maxsize = std::max( maxsize, totalKey.size() );
      // Compute the overall total
      size_t ipart = 0;
      float total = 0;
      float total2 = 0;
      for ( auto ip : m_partitionTimers )
      {
        total += ip.second;
        if ( ipart != 0 && ipart+1 != m_partitionTimers.size() ) total2 += ip.second;
        ipart++;
      }      
      // Dump individual partition timers and the overall total 
      for ( auto ip : m_partitionTimers )
        std::cout << std::setw(maxsize) << ip.first << " : " 
                  << std::fixed << std::setw(8) << ip.second << " sec" << std::endl;
      std::cout << std::setw(maxsize) << totalKey << " : " 
                << std::fixed << std::setw(8) << total << " sec" << std::endl;
      std::cout << std::setw(maxsize) << total2Key << " : " 
                << std::fixed << std::setw(8) << total2 << " sec" << std::endl;
    }

  private:

    Timer<TIMERTYPE> m_timer;
    std::string m_active;
    std::map< std::string, float > m_partitionTimers;
    std::map< std::string, uint32_t > m_partitionIds;

  };

}

#endif // MGONGPUTIMERMAP_H
