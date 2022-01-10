#ifndef EventStatistics_H
#define EventStatistics_H 1

//#include "mgOnGpuConfig.h"

#include <cmath>
//#include <sstream>

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
    int nevtOK; // number of events used, where ME is not nan (i.e. total - nan)
    int nevtNAN; // number of events used, where ME is nan
    double minME; // minimum matrix element
    double maxME; // maximum matrix element
    double minWG; // minimum sampling weight
    double maxWG; // maximum sampling weight
    double sumMEdiff; // sum of diff to mean for matrix element
    double sumWGdiff; // sum of diff to mean for sampling weight
    double sqsMEdiff; // squared sum of diff to mean for matrix element
    double sqsWGdiff; // squared sum of diff to mean for sampling weight
    double meanME() const { return minME + sumMEdiff / nevtOK; } // mean matrix element
    double meanWG() const { return minWG + sumWGdiff / nevtOK; } // mean sampling weight
    double stdME() const { return std::sqrt( sqsMEdiff / nevtOK ); } // standard deviation matrix element
    double stdWG() const { return std::sqrt( sqsWGdiff / nevtOK ); } // standard deviation sampling weight
  };

  //--------------------------------------------------------------------------

}

#endif // EventStatistics_H
