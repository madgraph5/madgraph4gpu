C     
C     If VECSIZE_MEMMAX is greater than 1, a vector API is used:
C     this is designed for offloading MEs to GPUs or vectorized C++,
C     but it can also be used for computing MEs in Fortran.
C     If VECSIZE_MEMMAX equals 1, the old scalar API is used:
C     this can only be used for computing MEs in Fortran.
C     
C     Fortran arrays in the vector API can hold up to VECSIZE_MEMMAX
C     events and are statically allocated at compile time.
C     The constant value of VECSIZE_MEMMAX is fixed at codegen time
C     (output madevent ... --vector_size=<VECSIZE_MEMMAX>).
C     
C     While the arrays can hold up to VECSIZE_MEMMAX events,
C     only VECSIZE_USED (<= VECSIZE_MEMAMX) are used in Fortran loops.
C     The value of VECSIZE_USED can be chosen at runtime
C     (typically 8k-16k for GPUs, 16-32 for vectorized C++).
C     
C     The value of VECSIZE_USED represents the number of events
C     handled by one call to the Fortran/cudacpp "bridge".
C     This is not necessarily the number of events which are
C     processed in lockstep within a single SIMD vector on CPUs
C     or within a single "warp" of threads on GPUs. These parameters
C     are internal to the cudacpp bridge and need not be exposed
C     to the Fortran program which calls the cudacpp bridge.
C     
C     NB: THIS FILE CANNOT CONTAIN #ifdef DIRECTIVES
C     BECAUSE IT DOES NOT GO THROUGH THE CPP PREPROCESSOR
C     (see https://github.com/madgraph5/madgraph4gpu/issues/458).
C     
      INTEGER WRAP_SIZE
      PARAMETER (WRAP_SIZE=32)
      INTEGER NB_WRAP
      PARAMETER (NB_WRAP=1)
      INTEGER VECSIZE_MEMMAX
      PARAMETER (VECSIZE_MEMMAX=32)
