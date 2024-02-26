
#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void ixxxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void ixxxxx(
#endif
               const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void ipzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void ipzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void imzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void imzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void ixzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void ixzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void vxxxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void vxxxxx( 
#endif
               const fptype_sv* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void sxxxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void sxxxxx( 
#endif
               const fptype_sv* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void oxxxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void oxxxxx( 
#endif
               const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void opzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void opzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void omzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void omzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC INLINE
  void oxzxxx( T_Acc const &acc,
#else
  __device__ INLINE
  void oxzxxx( 
#endif
               const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#if !defined(ALPAKA) && !defined(__CUDACC__)
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------
