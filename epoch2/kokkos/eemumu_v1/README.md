## Epoch2 Kokkos ee->mumu V1

This version was derived from the epoch1 Kokkos version. Key features:
* uses shared memory objects for momenta, weights, etc. 
* dedicated Random Number generator from Kokkos WITH memory cache
* new style rambo with `get_initial_momenta` and `get_final_momenta`, uses a `parallel_for(TeamPolicy)` inside
* sigmaKin uses monolithic `parallel_for(TeamPolicy)`

Build the source against CUDA and/or OpenMP by editing `src/Makefile` and `SubProcesses/P1_sigma_sm_epem_mupmum/Makefile`. Need to set `KOKKOSPATH_CUDA`, `KOKKOSPATH_OMP` and `CUDA_ARCH_NUM` (compute architecture).
Then run `make` inside `SubProcesses/P1_sigma_sm_epem_mupmum/Makefile`
One can also only build the CUDA or OpenMP versions with `make cuda` or `make openmp`.

Profiling this shows:
* `generate_random()` takes 0.077ms
* `get_final_momenta()` takes 3ms
* `sigmaKin()` takes 11.1ms
* Total time for 1 iteration: 15ms

## GTX 2070 running Kokkos CUDA Version

```
> ./check.exe -p 1000 256 100
***********************************
NumIterations         = 100
NumThreadsPerBlock    = 256
NumBlocksPerGrid      = 1000
-----------------------------------
TotalTimeInWaveFuncs  = 6.903039e-01 sec
MeanTimeInWaveFuncs   = 6.903039e-03 sec
StdDevTimeInWaveFuncs = 7.095544e-03 sec
MinTimeInWaveFuncs    = 5.635369e-03 sec
MaxTimeInWaveFuncs    = 7.667044e-02 sec
-----------------------------------
ProcessID:            = 13660
NProcesses            = 1
NumMatrixElements     = 25600000
MatrixElementsPerSec  = 3.708512e+07 sec^-1
***********************************
NumMatrixElements     = 25600000
MeanMatrixElemValue   = 1.252590e-02 GeV^0
StdErrMatrixElemValue = 1.724636e-06 GeV^0
StdDevMatrixElemValue = 8.726045e-03 GeV^0
MinMatrixElemValue    = 6.071582e-03 GeV^0
MaxMatrixElemValue    = 3.374925e-02 GeV^0
***********************************
fill_random_numbers   = 8.480381e-05 +/- 1.310079e-04 seconds
copy momenta          = 1.116470e-06 +/- 1.918433e-07 seconds
get final momenta     = 1.775829e-03 +/- 4.024518e-04 seconds
sigmaKin              = 6.903039e-03 +/- 7.095544e-03 seconds
copy momenta          = 1.374620e-06 +/- 1.585068e-07 seconds
copy matrix_element   = 5.935001e-07 +/- 2.729234e-08 seconds
full iteration        = 9.802165e-03 +/- 7.564926e-03 seconds
total time: 9.914289e-01  (compared to 1.6 seconds in epoch2/cuda/eemumu)
```

```
> ./check.exe -p 10000 256 100
***********************************
NumIterations         = 100
NumThreadsPerBlock    = 256
NumBlocksPerGrid      = 10000
-----------------------------------
TotalTimeInWaveFuncs  = 6.020865e+00 sec
MeanTimeInWaveFuncs   = 6.020865e-02 sec
StdDevTimeInWaveFuncs = 3.269770e-02 sec
MinTimeInWaveFuncs    = 5.622345e-02 sec
MaxTimeInWaveFuncs    = 3.854428e-01 sec
-----------------------------------
ProcessID:            = 13685
NProcesses            = 1
NumMatrixElements     = 256000000
MatrixElementsPerSec  = 4.251881e+07 sec^-1
***********************************
NumMatrixElements     = 256000000
MeanMatrixElemValue   = 4.002410e-03 GeV^0
StdErrMatrixElemValue = 6.227252e-07 GeV^0
StdDevMatrixElemValue = 9.963603e-03 GeV^0
MinMatrixElemValue    = 6.071582e-03 GeV^0
MaxMatrixElemValue    = 3.374926e-02 GeV^0
***********************************
fill_random_numbers   = 8.707275e-04 +/- 7.314602e-04 seconds
copy momenta          = 1.546860e-06 +/- 2.731407e-07 seconds
get final momenta     = 1.366333e-02 +/- 1.380719e-03 seconds
sigmaKin              = 6.020865e-02 +/- 3.269770e-02 seconds
copy momenta          = 1.912520e-06 +/- 5.860882e-07 seconds
copy matrix_element   = 6.390500e-07 +/- 3.886275e-08 seconds
full iteration        = 8.504003e-02 +/- 3.492964e-02 seconds
total time: 8.556520e+00
```

## GTX 2070 running OpenMP version
```
-> export OMP_NUM_THREADS=64;./ocheck.exe -p $(( 256000 / $OMP_NUM_THREADS )) $OMP_NUM_THREADS 100
Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false
Kokkos::OpenMP::initialize WARNING: You are likely oversubscribing your CPU cores.
                                    Detected: 12 cores per node.
                                    Detected: 1 MPI_ranks per node.
                                    Requested: 64 threads per process.
Opened slha file ../../Cards/param_card.dat for reading
***********************************
NumIterations         = 100
NumThreadsPerBlock    = 64
NumBlocksPerGrid      = 4000
-----------------------------------
TotalTimeInWaveFuncs  = 5.726495e+01 sec
MeanTimeInWaveFuncs   = 5.726495e-01 sec
StdDevTimeInWaveFuncs = 8.226839e-02 sec
MinTimeInWaveFuncs    = 4.799600e-01 sec
MaxTimeInWaveFuncs    = 1.245278e+00 sec
-----------------------------------
ProcessID:            = 34985
NProcesses            = 1
NumMatrixElements     = 25600000
MatrixElementsPerSec  = 4.470448e+05 sec^-1
***********************************
NumMatrixElements     = 25600000
MeanMatrixElemValue   = 1.252550e-02 GeV^0
StdErrMatrixElemValue = 1.725479e-06 GeV^0
StdDevMatrixElemValue = 8.730309e-03 GeV^0
MinMatrixElemValue    = 6.071582e-03 GeV^0
MaxMatrixElemValue    = 3.374926e-02 GeV^0
***********************************
fill_random_numbers   = 5.399807e-01 +/- 4.495470e-02 seconds
get_initial_momenta   = 5.360466e-01 +/- 4.466969e-02 seconds
get_final_momenta     = 5.406816e-01 +/- 3.863776e-02 seconds
copy weights          = 8.439896e-07 +/- 1.146274e-07 seconds
copy momenta          = 5.732799e-07 +/- 9.457919e-08 seconds
sigmaKin              = 5.726495e-01 +/- 8.226839e-02 seconds
copy matrix_element   = 1.155050e-06 +/- 1.516596e-07 seconds
full iteration        = 2.190087e+00 +/- 1.228588e-01 seconds
total time: 2.190111e+02
```


## GTX 2070 running CUDA version in `epoch2/cuda/eemumu`
```
-> ./gcheck.exe -p 10000 256 100
***********************************************************************
NumBlocksPerGrid           = 10000
NumThreadsPerBlock         = 256
NumIterations              = 100
-----------------------------------------------------------------------
FP precision               = DOUBLE (nan=0)
Complex type               = THRUST::COMPLEX
RanNumb memory layout      = AOSOA[4]
Momenta memory layout      = AOSOA[4]
Wavefunction GPU memory    = LOCAL
Random number generation   = CURAND DEVICE (CUDA code)
-----------------------------------------------------------------------
NumberOfEntries            = 100
TotalTime[Rnd+Rmb+ME] (123)= ( 1.093454e+01                 )  sec
TotalTime[Rambo+ME]    (23)= ( 1.089835e+01                 )  sec
TotalTime[RndNumGen]    (1)= ( 3.618864e-02                 )  sec
TotalTime[Rambo]        (2)= ( 6.780948e+00                 )  sec
TotalTime[MatrixElems]  (3)= ( 4.117404e+00                 )  sec
MeanTimeInMatrixElems      = ( 4.117404e-02                 )  sec
[Min,Max]TimeInMatrixElems = [ 4.106871e-02 ,  4.210467e-02 ]  sec
-----------------------------------------------------------------------
TotalEventsComputed        = 256000000
EvtsPerSec[Rnd+Rmb+ME](123)= ( 2.341205e+07                 )  sec^-1
EvtsPerSec[Rmb+ME]     (23)= ( 2.348979e+07                 )  sec^-1
EvtsPerSec[MatrixElems] (3)= ( 6.217510e+07                 )  sec^-1
***********************************************************************
NumMatrixElements(notNan)  = 256000000
MeanMatrixElemValue        = ( 1.371753e-02 +- 5.125405e-07 )  GeV^0
[Min,Max]MatrixElemValue   = [ 6.071581e-03 ,  3.374941e-02 ]  GeV^0
StdDevMatrixElemValue      = ( 8.200648e-03                 )  GeV^0
MeanWeight                 = ( 4.515827e-01 +- 0.000000e+00 )
[Min,Max]Weight            = [ 4.515827e-01 ,  4.515827e-01 ]
StdDevWeight               = ( 0.000000e+00                 )
***********************************************************************
00 CudaFree :     0.101069 sec
0a ProcInit :     0.000142 sec
0b MemAlloc :     0.070105 sec
0c GenCreat :     0.007103 sec
0d SGoodHel :     0.245752 sec
1a GenSeed  :     0.000105 sec
1b GenRnGen :     0.036084 sec
2a RamboIni :     0.000382 sec
2b RamboFin :     0.000195 sec
2c CpDTHwgt :     1.774336 sec
2d CpDTHmom :     5.006034 sec
3a SigmaKin :     3.803070 sec
3b CpDTHmes :     0.314333 sec
4a DumpLoop :     1.314959 sec
8a CompStat :     1.146471 sec
9a GenDestr :     0.000025 sec
9b DumpScrn :     0.000081 sec
9c DumpJson :     0.000001 sec
TOTAL       :    13.820248 sec
TOTAL (123) :    10.934539 sec
TOTAL  (23) :    10.898350 sec
TOTAL   (1) :     0.036189 sec
TOTAL   (2) :     6.780947 sec
TOTAL   (3) :     4.117404 sec
***********************************************************************
```
