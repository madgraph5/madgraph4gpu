C
C All Fortran arrays for cudacpp offload can hold up to NB_PAGE_MAX events 
C These arrays are allocated at compile time (NB_PAGE_MAX is a constant)
C      
C Of these, only NB_PAGE_LOOP events are used in Fortran loops (#455)
C The value of NB_PAGE_LOOP can be chosen at runtime
C (very large for CUDA/GPUs, smaller for C++/CPUs using SIMD)
C
C NB: THIS FILE CANNOT CONTAIN #ifdef DIRECTIVES BECAUSE IT DOES NOT GO THROUGH THE CPP PREPROCESSOR (#458)
C 
C
      INTEGER NB_PAGE_MAX
      PARAMETER (NB_PAGE_MAX=16384) ! NB: 16k events per GPU grid is the minimum required to fill a V100 GPU
c     PARAMETER (NB_PAGE_MAX=32) ! NB: workaround for out-of-memory on Juwels: 32 is enough for no-CUDA builds (issue #498)
      INTEGER NB_PAGE_LOOP
      COMMON/TO_NB_PAGE/NB_PAGE_LOOP
