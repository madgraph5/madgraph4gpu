C
C All Fortran arrays for cudacpp offload can hold up to VECSIZE_MEMMAX events 
C These arrays are allocated at compile time (VECSIZE_MEMMAX is a constant)
C      
C Of these, only VECSIZE_USED events are used in Fortran loops (#455)
C The value of VECSIZE_USED can be chosen at runtime
C (very large for CUDA/GPUs, smaller for C++/CPUs using SIMD)
C
C NB: THIS FILE CANNOT CONTAIN #ifdef DIRECTIVES BECAUSE IT DOES NOT GO THROUGH THE CPP PREPROCESSOR (#458)
C 
C
      INTEGER VECSIZE_MEMMAX
      PARAMETER (VECSIZE_MEMMAX=16384) ! NB: 16k events per GPU grid is the minimum required to fill a V100 GPU
c     PARAMETER (VECSIZE_MEMMAX=32) ! NB: workaround for out-of-memory on Juwels: 32 is enough for no-CUDA builds (issue #498)
      INTEGER VECSIZE_USED
      COMMON/TO_VECSIZE/VECSIZE_USED
