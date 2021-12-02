
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
  void ixxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void ipzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void imzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__ INLINE
  void ixzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void vxxxxx( const fptype* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void sxxxxx( const fptype* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__ INLINE
  void opzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__ INLINE
  void omzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__ INLINE
  void oxzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------
