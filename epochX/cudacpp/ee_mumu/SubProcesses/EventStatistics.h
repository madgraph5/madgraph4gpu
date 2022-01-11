#ifndef EventStatistics_H
#define EventStatistics_H 1

#include <algorithm>
#include <cmath>
#include <limits>

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  // The EventStatistics struct is used to accumulate running aggregates of event statistics.
  // This include the process cross section and the process maximum weight.
  // One important case of EventStatistics is the "gridpack" result set, which is
  // the output of the "integration" step and the input to "unweighted event generation" step.
  // In first approximation, the process cross section and maximum weight are just the mean ME and maximum ME;
  // eventually, sampling weights (e.g. from Rambo) must also be taken into account in this calculation.
  struct EventStatistics
  {
    size_t nevtALL; // total number of events used
    size_t nevtABN; // number of events used, where ME is abnormal (nevtABN <= nevtALL)
    size_t nevtZERO; // number of not-abnormal events used, where ME is zero (nevtZERO <= nevtOK)
    double minME; // minimum matrix element
    double maxME; // maximum matrix element
    double minWG; // minimum sampling weight
    double maxWG; // maximum sampling weight
    double sumMEdiff; // sum of diff to mean for matrix element
    double sumWGdiff; // sum of diff to mean for sampling weight
    double sqsMEdiff; // squared sum of diff to mean for matrix element
    double sqsWGdiff; // squared sum of diff to mean for sampling weight
    size_t nevtOK() const { return nevtALL - nevtABN; } // number of events used, where ME is not abnormal
    double meanME() const { return minME + sumMEdiff / nevtOK(); } // mean matrix element
    double meanWG() const { return minWG + sumWGdiff / nevtOK(); } // mean sampling weight
    double stdME() const { return std::sqrt( sqsMEdiff / nevtOK() ); } // standard deviation matrix element
    double stdWG() const { return std::sqrt( sqsWGdiff / nevtOK() ); } // standard deviation sampling weight
    // Constructor
    EventStatistics()
      : nevtALL( 0 )
      , nevtABN( 0 )
      , nevtZERO( 0 )
      , minME( std::numeric_limits<double>::max() )
      , maxME( std::numeric_limits<double>::min() )
      , minWG( std::numeric_limits<double>::max() )
      , maxWG( std::numeric_limits<double>::min() )
      , sumMEdiff( 0 )
      , sumWGdiff( 0 )
      , sqsMEdiff( 0 )
      , sqsWGdiff( 0 ){}
    // Combine two EventStatistics
    EventStatistics& operator+=( const EventStatistics& stats )
    {
      const EventStatistics s1 = *this; // temporary const copy
      const EventStatistics& s2 = stats;
      EventStatistics& sum = *this;
      sum.nevtALL = s1.nevtALL + s2.nevtALL;
      sum.nevtABN = s1.nevtABN + s2.nevtABN;
      sum.nevtZERO = s1.nevtZERO + s2.nevtZERO;
      sum.minME = std::min( s1.minME, s2.minME );
      sum.maxME = std::max( s1.maxME, s2.maxME );
      sum.minWG = std::min( s1.minWG, s2.minWG );
      sum.maxWG = std::max( s1.maxWG, s2.maxWG );
      sum.sumMEdiff =
        s1.sumMEdiff + s1.nevtOK() * ( s1.minME - sum.minME ) +
        s2.sumMEdiff + s2.nevtOK() * ( s2.minME - sum.minME );
      sum.sumWGdiff =
        s1.sumWGdiff + s1.nevtOK() * ( s1.minWG - sum.minWG ) +
        s2.sumWGdiff + s2.nevtOK() * ( s2.minWG - sum.minWG );
      sum.sqsMEdiff =
        s1.sqsMEdiff + s1.nevtOK() * std::pow( s1.meanME() - sum.minME, 2 ) +
        s2.sqsMEdiff + s2.nevtOK() * std::pow( s2.meanME() - sum.minME, 2 );
      sum.sqsWGdiff =
        s1.sqsWGdiff + s1.nevtOK() * std::pow( s1.meanWG() - sum.minWG, 2 ) +
        s2.sqsWGdiff + s2.nevtOK() * std::pow( s2.meanWG() - sum.minWG, 2 );
      return sum;
    }
  };

  //--------------------------------------------------------------------------

}

#endif // EventStatistics_H
