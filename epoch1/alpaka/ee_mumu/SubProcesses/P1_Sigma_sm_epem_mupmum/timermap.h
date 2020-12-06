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
      const std::string totalKey = "TOTAL      "; // "TOTAL(ANY) "?
      //const std::string totalBut2Key = "TOTAL(n-2) ";
      const std::string total123Key = "TOTAL(123) ";
      const std::string total23Key = "TOTAL(23)  ";
      const std::string total3Key = "TOTAL(3)   ";
      size_t maxsize = 0;
      for ( auto ip : m_partitionTimers )
        maxsize = std::max( maxsize, ip.first.size() );
      maxsize = std::max( maxsize, totalKey.size() );
      // Compute the overall total
      size_t ipart = 0;
      float total = 0;
      //float totalBut2 = 0;
      float total123 = 0;
      float total23 = 0;
      float total3 = 0;
      for ( auto ip : m_partitionTimers )
      {
        total += ip.second;
        //if ( ipart != 0 && ipart+1 != m_partitionTimers.size() ) totalBut2 += ip.second;
        if ( ip.first[0] == '1' || ip.first[0] == '2' || ip.first[0] == '3' ) total123 += ip.second;
        if ( ip.first[0] == '2' || ip.first[0] == '3' ) total23 += ip.second;
        if ( ip.first[0] == '3' ) total3 += ip.second;
        ipart++;
      }      
      // Dump individual partition timers and the overall total 
      for ( auto ip : m_partitionTimers )
        std::cout << std::setw(maxsize) << ip.first << " : " 
                  << std::fixed << std::setw(8) << ip.second << " sec" << std::endl;
      std::cout << std::setw(maxsize) << totalKey << " : " 
                << std::fixed << std::setw(8) << total << " sec" << std::endl;
      //std::cout << std::setw(maxsize) << totalBut2Key << " : " 
      //          << std::fixed << std::setw(8) << totalBut2 << " sec" << std::endl;
      std::cout << std::setw(maxsize) << total123Key << " : " 
                << std::fixed << std::setw(8) << total123 << " sec" << std::endl;
      std::cout << std::setw(maxsize) << total23Key << " : " 
                << std::fixed << std::setw(8) << total23 << " sec" << std::endl;
      std::cout << std::setw(maxsize) << total3Key << " : " 
                << std::fixed << std::setw(8) << total3 << " sec" << std::endl;
    }

  private:

    Timer<TIMERTYPE> m_timer;
    std::string m_active;
    std::map< std::string, float > m_partitionTimers;
    std::map< std::string, uint32_t > m_partitionIds;

  };

}

#endif // MGONGPUTIMERMAP_H
