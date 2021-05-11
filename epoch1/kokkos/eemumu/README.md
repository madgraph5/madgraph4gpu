## Epoch1 Kokkos ee->mumu 

This version was derived from the epoch0 CUDA version. Key features:
* dedicated Random Number generator from Kokkos, memory generated each time
* old style rambo, which creates the momenta memory, uses a `parallel_for(RangePolicy)` inside
* sigmaKin uses `parallel_for(RangePolicy)` as well,
* not much shared memory use, frequent creation/destruction

Profiling this shows:
* `generate_random()` takes 0.5ms
* `sigmaKin()` takes 20.4ms
* Total time for 1 iteration: 29ms


## GTX 2070 

```
> ./check.exe -p 1000 256 100
Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
  In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
  For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
  For unit testing set OMP_PROC_BIND=false
Kokkos::Cuda::initialize WARNING: Cuda is allocating into UVMSpace by default
                                  without setting CUDA_LAUNCH_BLOCKING=1.
                                  The code must call Cuda().fence() after each kernel
                                  or will likely crash when accessing data on the host.
Kokkos::Cuda::initialize WARNING: Cuda is allocating into UVMSpace by default
                                  without setting CUDA_MANAGED_FORCE_DEVICE_ALLOC=1 or 
                                  setting CUDA_VISIBLE_DEVICES.
                                  This could on multi GPU systems lead to severe performance
                                  penalties.
Opened slha file ../../Cards/param_card.dat for reading
***********************************
NumIterations         = 100
NumThreadsPerBlock    = 256
NumBlocksPerGrid      = 1000
-----------------------------------
NumberOfEntries       = 100
TotalTimeInWaveFuncs  = 2.031924e+00 sec
MeanTimeInWaveFuncs   = 2.031924e-02 sec
StdDevTimeInWaveFuncs = 4.430294e-03 sec
MinTimeInWaveFuncs    = 1.929230e-02 sec
MaxTimeInWaveFuncs    = 4.852227e-02 sec
-----------------------------------
ProcessID:            = 36035
NProcesses            = 1
NumMatrixElements     = 25600000
MatrixElementsPerSec  = 1.259890e+07 sec^-1
***********************************
NumMatrixElements     = 25600000
MeanMatrixElemValue   = 1.373249e-02 GeV^0
StdErrMatrixElemValue = 1.622932e-06 GeV^0
StdDevMatrixElemValue = 8.211460e-03 GeV^0
MinMatrixElemValue    = 6.071582e-03 GeV^0
MaxMatrixElemValue    = 3.374915e-02 GeV^0
total time: 2.909314e+00 (compared to 1.8 seconds in epoch2/cuda/eemumu)
```
