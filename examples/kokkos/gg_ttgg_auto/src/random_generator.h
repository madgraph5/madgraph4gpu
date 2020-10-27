#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

// Inspired by this example:
// https://github.com/kokkos/kokkos/blob/master/example/tutorial/Algorithms/01_random_numbers/random_numbers.cpp

// A Functor for generating uint64_t random numbers templated on the
// GeneratorPool type
template <class GeneratorPool,typename ExecSpace>
struct generate_random {
  // Output View for the random numbers
  Kokkos::View<double**,ExecSpace> vals;

  // The GeneratorPool
  GeneratorPool rand_pool;

  int samples;

  // Initialize all members
  generate_random(Kokkos::View<double**,ExecSpace> vals_, GeneratorPool rand_pool_,
                  int samples_)
      : vals(vals_), rand_pool(rand_pool_), samples(samples_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    // Get a random number state from the pool for the active thread
    typename GeneratorPool::generator_type rand_gen = rand_pool.get_state();

    // Draw samples numbers from the pool as urand64 between 0 and
    // rand_pool.MAX_URAND64 Note there are function calls to get other type of
    // scalars, and also to specify Ranges or get a normal distributed float.
    for (int k = 0; k < samples; k++)
      vals(i,k) = rand_gen.frand();

    // Give the state back, which will allow another thread to acquire it
    rand_pool.free_state(rand_gen);
  }
};


template <typename ExecSpace>
Kokkos::View<double**,ExecSpace> get_random_numbers(const int size,const int samples,const uint64_t random_seed=5374857){
  Kokkos::View<double**,ExecSpace> rns("rns",size,samples);

  Kokkos::Random_XorShift64_Pool<> rand_pool64(random_seed);

  Kokkos::parallel_for("generate_random_numbers",Kokkos::RangePolicy<ExecSpace>(0,size),
    generate_random<Kokkos::Random_XorShift64_Pool<>,ExecSpace>(rns, rand_pool64, samples));

  return rns;
}



#endif