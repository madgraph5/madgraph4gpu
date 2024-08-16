// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Aug 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#ifndef MGONGPUTIMERS_H
#define MGONGPUTIMERS_H 1

#include <cassert>
#include <chrono>
#include <iostream>
#include <type_traits>

namespace mgOnGpu
{

  // ---------------------------------------------------------------------------
  
  // Default ("old") timers based on std::chrono clocks
  // Based on the previous timer.h header by S.Roiser, O. Mattelaer and A. Valassi
  // Template argument T can be any of high_resolution_clock, steady_clock, system_clock
  // See https://www.modernescpp.com/index.php/the-three-clocks
  // See https://codereview.stackexchange.com/questions/196245/extremely-simple-timer-class-in-c  
  template<typename T>
  class ChronoTimer
  {
  public:
    ChronoTimer();
    virtual ~ChronoTimer() {}
    void start();
    void stop();
    float getDurationSeconds();
  private:
    std::chrono::duration<float> m_duration;
    bool m_started;
    typedef typename T::time_point TTP;
    TTP m_startTime;
  };

  template<typename T>
  ChronoTimer<T>::ChronoTimer()
    : m_duration( 0 )
    , m_started( false )
    , m_startTime()
  {
    static_assert( std::is_same<T, std::chrono::high_resolution_clock>::value ||
                   std::is_same<T, std::chrono::steady_clock>::value ||
                   std::is_same<T, std::chrono::system_clock>::value );
  }

  template<typename T>
  void
  ChronoTimer<T>::start()
  {
    assert( !m_started );
    m_started = true;
    m_startTime = T::now();
  }

  template<typename T>
  void
  ChronoTimer<T>::stop()
  {
    assert( m_started );
    m_started = false;
    m_duration += T::now() - m_startTime;
  }

  template<typename T>
  float
  ChronoTimer<T>::getDurationSeconds()
  {
    assert( !m_started );
    return m_duration.count();
  }

  // ---------------------------------------------------------------------------

}
#endif // MGONGPUTIMERS_H
