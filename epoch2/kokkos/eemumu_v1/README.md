## Epoch2 Kokkos ee->mumu V1

This version was derived from the epoch1 Kokkos version. Key features:
* uses shared memory objects for momenta, weights, etc. 
* dedicated Random Number generator from Kokkos WITH memory cache
* new style rambo with `get_initial_momenta` and `get_final_momenta`, momenta stored in cache, uses a `parallel_for(TeamPolicy)` inside
* sigmaKin uses monolithic `parallel_for(TeamPolicy)`

Profiling this shows:
* `get_momenta()` takes 1.5ms
* `sigmaKin()` takes 20.4ms
* Total time for 1 iteration: 29.5ms
