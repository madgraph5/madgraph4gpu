C     
C     Offloading ME to GPUs or vectorized C++ needs a vector API.
C     Fortran arrays in this API can hold up to VECSIZE_MEMMAX events.
C     These arrays are statically allocated at compile time.
C     The constant value of VECSIZE_MEMMAX is fixed at codegen time
C     (output madevent ... --vector_size=<VECSIZE_MEMMAX>).
C     
C     While the arrays can hold up to VECSIZE_MEMMAX events,
C     only VECSIZE_USED (<= VECSIZE_MEMAMX) are used in Fortran loops.
C     The value of VECSIZE_USED can be chosen at runtime
C     (typically 8k-16k for GPUs, 16-32 for vectorized C++).
C     
C     NB: THIS FILE CANNOT CONTAIN #ifdef DIRECTIVES
C     BECAUSE IT DOES NOT GO THROUGH THE CPP PREPROCESSOR
C     (see https://github.com/madgraph5/madgraph4gpu/issues/458).
C     
      INTEGER VECSIZE_MEMMAX
      PARAMETER (VECSIZE_MEMMAX=16384) ! NB: 16k events per GPU grid is the minimum required to fill a V100 GPU
c     PARAMETER (VECSIZE_MEMMAX=32) ! NB: workaround for out-of-memory on Juwels: 32 is enough for no-CUDA builds (issue #498)
