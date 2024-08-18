// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
//==========================================================================
// Created by: S. Roiser (Feb 2020) for the MG5aMC CUDACPP plugin [old chrono timer, old API].
// Further modified by: O. Mattelaer, S. Roiser, A. Valassi (2020-2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Created by: A. Valassi (Aug 2024) for the MG5aMC CUDACPP plugin [new chrono timer, new API, add rdtsc timer].
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.
//==========================================================================

#ifndef MGONGPUTIMER_H
#define MGONGPUTIMER_H 1

#include <cassert>
#include <chrono>
#include <iostream>
#include <ratio>
#include <type_traits>

namespace mgOnGpu
{

  // ---------------------------------------------------------------------------
  
  // ChronoTimer: default ("old") timers based on std::chrono clocks
  // With respect to the original Timer class, this uses a new implementation with nanosecond counts
  // With respect to the original Timer class, this also uses a new API with explicit start/stop
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
    float getDurationSeconds( bool allowRunning = false ); // by default, assert that the timer is not running
  private:
    typedef std::nano RATIO;
    typedef std::chrono::duration<uint64_t, RATIO> DURATION;
    typedef std::chrono::time_point<T, DURATION> TIMEPOINT;
    DURATION m_duration;
    bool m_started;
    TIMEPOINT m_startTime;
  };

  template<typename T>
  inline
  ChronoTimer<T>::ChronoTimer()
    : m_duration()
    , m_started( false )
    , m_startTime()
  {
    static_assert( std::is_same<T, std::chrono::high_resolution_clock>::value ||
                   std::is_same<T, std::chrono::steady_clock>::value ||
                   std::is_same<T, std::chrono::system_clock>::value );
  }

  template<typename T>
  inline
  void
  ChronoTimer<T>::start()
  {
    assert( !m_started );
    m_started = true;
    m_startTime = T::now();
  }

  template<typename T>
  inline
  void
  ChronoTimer<T>::stop()
  {
    assert( m_started );
    m_started = false;
    m_duration += T::now() - m_startTime;
  }

  template<typename T>
  inline
  float
  ChronoTimer<T>::getDurationSeconds( bool allowRunning )
  {
    if( allowRunning ) stop(); // (old timer behaviour) compute m_duration and allow next start() call
    assert( !m_started );
    auto count = m_duration.count();
    if( allowRunning ) m_duration = DURATION::zero(); // (old timer behaviour) reset m_duration
    return count * (float)RATIO::num / RATIO::den;
  }

  // ---------------------------------------------------------------------------
  
  // RdtscTimer: faster ("new") *EXPERIMENTAL* timers based on rdtsc
  // The rdtsc() call is derived from the TSCNS class (https://github.com/MengRao/tscns)
  // The conversion of rdtsc counts to seconds is calibrated on the average frequency during the timer lifetime
  // See https://stackoverflow.com/q/76063685 and the Intel 64 and IA-32 Architectures Software Developer’s Manual
  // (https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html, June 2024):
  // "To determine average processor clock frequency, Intel recommends the use of performance monitoring
  // logic to count processor core clocks over the period of time for which the average is required."
  class RdtscTimer
  {
  public:
    RdtscTimer();
    virtual ~RdtscTimer() {}
    void start();
    void stop();
    float getDurationSeconds( bool allowRunning = false ); // by default, assert that the timer is not running
    static uint64_t rdtsc();
  private:
    uint64_t m_duration;
    bool m_started;
    uint64_t m_startCount;
    ChronoTimer<std::chrono::high_resolution_clock> m_ctorTimer;
    uint64_t m_ctorCount;
  };

  inline
  uint64_t
  RdtscTimer::rdtsc()
  {
#if defined( __x86_64__ )
    return __builtin_ia32_rdtsc();
#else
#error "rdtsc is not defined for this platform yet"
#endif
  }

  inline
  void
  RdtscTimer::start()
  {
    assert( !m_started );
    m_started = true;
    m_startCount = rdtsc();
  }

  inline
  RdtscTimer::RdtscTimer()
    : m_duration( 0 )
    , m_started( false )
    , m_startCount( 0 )
    , m_ctorTimer()
    , m_ctorCount( 0 )
  {
    m_ctorTimer.start();
    m_ctorCount = rdtsc();
  }

  inline
  void
  RdtscTimer::stop()
  {
    assert( m_started );
    m_started = false;
    m_duration += rdtsc() - m_startCount;
  }

  inline
  float
  RdtscTimer::getDurationSeconds( bool allowRunning )
  {
    if( allowRunning ) stop(); // (old timer behaviour) compute m_duration and allow next start() call
    assert( !m_started );
    m_ctorTimer.stop();
    float secPerCount = m_ctorTimer.getDurationSeconds() / ( rdtsc() - m_ctorCount );
    m_ctorTimer.start(); // just in case getDurationSeconds() is called again... (e.g. if allowRunning is true)
    auto count = m_duration;
    if( allowRunning ) m_duration = 0; // (old timer behaviour) reset m_duration
    return count * secPerCount;
  }

  // ---------------------------------------------------------------------------

}
#endif // MGONGPUTIMER_H
