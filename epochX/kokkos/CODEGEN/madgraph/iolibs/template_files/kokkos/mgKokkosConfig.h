#ifndef MGKOKKOSCONFIG_H
#define MGKOKKOSCONFIG_H 1

namespace mgKokkos
{
  // --- Physics process-specific constants that are best declared at compile time

  constexpr int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
  constexpr int nw6 = %(wavefuncsize)d; // dimension of each wavefunction (see KEK 91-11)

  constexpr int npari = %(nincoming)d; // #particles in the initial state (incoming): e+ e-
  constexpr int nparf = %(noutcoming)d; // #particles in the final state (outgoing): mu+ mu-
  constexpr int npar = npari + nparf; // #particles in total (external): e+ e- -> mu+ mu-

  constexpr int nwf = %(nwavefunc)d; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)

  constexpr int ncomb = %(nbhel)d; // #helicity combinations: 16=2(spin up/down for fermions)**4(npar)

  // --- Platform-specific software implementation details

  // Maximum number of threads per block
  constexpr int ntpbMAX = 256;

}

#endif // MGKOKKOSCONFIG_H
