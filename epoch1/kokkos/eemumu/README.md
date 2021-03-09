## Epoch1 Kokkos ee->mumu 

This version was derived from the epoch0 CUDA version. Key features:
* old style rambo, which creates the momenta memory, uses a `parallel_for(RangePolicy)` inside
* sigmaKin uses `parallel_for(RangePolicy)` as well,
* not much shared memory use, frequent creation/destruction

Profiling this shows:
* `get_momenta()` takes 1.5ms
* `sigmaKin()` takes 20.4ms
* Total time for 1 iteration: 29ms
