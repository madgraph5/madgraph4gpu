## Epoch2 Kokkos ee->mumu V1

This version was derived from the epoch1 Kokkos version. Key features:
* uses shared memory objects for momenta, weights, etc. 
* dedicated Random Number generator from Kokkos WITH memory cache
* new style rambo with `get_initial_momenta` and `get_final_momenta`, momenta stored in cache, uses a `parallel_for(TeamPolicy)` inside
* sigmaKin uses monolithic `parallel_for(TeamPolicy)`

Profiling this shows:
* `generate_random()` takes 0.5ms
* `get_final_momenta()` takes 1.5ms
* `sigmaKin()` takes 19.6ms
* Total time for 1 iteration: 24ms

## GTX 2070 

```
> ./check.exe -p 1000 256 100
***********************************
NumIterations         = 100
NumThreadsPerBlock    = 256
NumBlocksPerGrid      = 1000
-----------------------------------
TotalTimeInWaveFuncs  = 2.020332e+00 sec
MeanTimeInWaveFuncs   = 2.020332e-02 sec
StdDevTimeInWaveFuncs = 8.694585e-04 sec
MinTimeInWaveFuncs    = 1.961946e-02 sec
MaxTimeInWaveFuncs    = 2.342804e-02 sec
-----------------------------------
ProcessID:            = 35882
NProcesses            = 1
NumMatrixElements     = 25600000
MatrixElementsPerSec  = 1.267118e+07 sec^-1
***********************************
NumMatrixElements     = 25600000
MeanMatrixElemValue   = 1.252599e-02 GeV^0
StdErrMatrixElemValue = 1.724525e-06 GeV^0
StdDevMatrixElemValue = 8.725484e-03 GeV^0
MinMatrixElemValue    = 6.071582e-03 GeV^0
MaxMatrixElemValue    = 3.374925e-02 GeV^0
***********************************
fill_random_numbers   = 1.081320e-04 +/- 1.316756e-04 seconds
copy momenta          = 1.222050e-06 +/- 2.832321e-07 seconds
get final momenta     = 1.624340e-03 +/- 2.191231e-05 seconds
sigmaKin              = 2.020332e-02 +/- 8.694585e-04 seconds
copy momenta          = 2.265440e-06 +/- 9.532062e-07 seconds
copy matrix_element   = 6.360600e-07 +/- 7.100688e-08 seconds
total time: 2.305009e+00  (compared to 1.8 seconds in epoch2/cuda/eemumu)
```
