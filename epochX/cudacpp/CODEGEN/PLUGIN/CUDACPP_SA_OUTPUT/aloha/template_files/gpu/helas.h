
  // =============================================================================
  // *** Generic pattern: kernelAccessFunction( buffer, additional_indexes ) ***
  // =============================================================================

  // Kernel access function (WITHOUT an explicit event number) for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the 4-momenta for one event or one SIMD vector of events
  // (Non-const memory access)
  __device__ inline
  fptype_sv& kernelAccessMomenta( fptype_sv* buffer,
                                  const int ip4
#ifdef __CUDACC__
                                  , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                  );

  // (Const memory access)
  __device__ inline
  const fptype_sv& kernelAccessConstMomenta( const fptype_sv* buffer,
                                             const int ip4
#ifdef __CUDACC__
                                             , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                             );

  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void ixxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void ipzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void imzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__ INLINE
  void ixzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void vxxxxx( const fptype_sv* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void sxxxxx( const fptype_sv* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxxxxx( const fptype_sv* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void opzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void omzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxzxxx( const fptype_sv* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
#ifdef __CUDACC__
               , const int ipar                // input: particle# out of npar
#endif
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------
