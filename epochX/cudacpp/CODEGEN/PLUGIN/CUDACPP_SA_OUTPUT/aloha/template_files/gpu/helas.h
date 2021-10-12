
#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ipzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void imzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void vxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv* vc,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void sxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype,                   // WARNING: "smass" unused (missing in Fortran)
               const int,                      // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[3],                // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxxxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void opzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void omzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxzxxx( const fptype_sv* allmomenta,    // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------
