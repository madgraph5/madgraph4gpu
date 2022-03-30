
  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------
  // Memory access functions
  SYCL_EXTERNAL
  inline const fptype& kernelAccessIp4IparIevt(
          const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
          const int ip4,
          const int ipar,
          const size_t ievt
          ) {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    const int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    return momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void ixxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL INLINE
  void ipzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL INLINE
  void imzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL INLINE
  void ixzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void vxxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void sxxxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype,                 // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               //const int,                    // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void oxxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL INLINE
  void opzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL INLINE
  void omzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void oxzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar                  // input: particle# out of npar
               ) ALWAYS_INLINE;

  //==========================================================================

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void ixxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fi = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    if ( fmass != 0. )
    {
      const fptype_sv pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
      if ( pp == 0. )
      {
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use sycl::abs!
        fptype sqm[2] = { fpsqrt( sycl::fabs( fmass ) ), 0. }; // possibility of negative fermion masses
        //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const int ip = ( 1 + nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        const int im = ( 1 - nh ) / 2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        fi[2] = cxmake( ip * sqm[ip], 0 );
        fi[3] = cxmake( im * nsf * sqm[ip], 0 );
        fi[4] = cxmake( ip * nsf * sqm[im], 0 );
        fi[5] = cxmake( im * sqm[im], 0 );
      }
      else
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype omega[2] = { fpsqrt( pvec0 + pp ), 0. };
        omega[1] = fmass / omega[0];
        const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype pp3 = fpmax( pp + pvec3, 0. );
        const cxtype chi[2] = { cxmake( fpsqrt ( pp3 * (fptype)0.5 / pp ), 0. ),
                                ( pp3 == 0. ?
                                  cxmake( -nh, 0. ) :
                                  cxmake( nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        fi[2] = sfomega[0] * chi[im];
        fi[3] = sfomega[0] * chi[ip];
        fi[4] = sfomega[1] * chi[im];
        fi[5] = sfomega[1] * chi[ip];
      }
    }
    else
    {
      const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. and pvec2 == 0. and pvec3 < 0. ),
                                          fptype_sv{0}, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
      const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ), cxternary( ( sqp0p3 == 0. ),
                                                                  cxmake( -(fptype)nhel * fpsqrt( 2. * pvec0 ), 0. ),
                                                                  cxmake( (fptype)nh * pvec1, pvec2 ) / sqp0p3 ) };
      if ( nh == 1 )
      {
        fi[2] = cxzero_sv();
        fi[3] = cxzero_sv();
        fi[4] = chi[0];
        fi[5] = chi[1];
      }
      else
      {
        fi[2] = chi[1];
        fi[3] = chi[0];
        fi[4] = cxzero_sv();
        fi[5] = cxzero_sv();
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL
  void ipzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fi = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fi[0] = cxmake( -pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv sqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
    fi[2] = fi[1];
    if( nh == 1 )
    {
      fi[3] = fi[1];
      fi[4] = sqp0p3;
    }
    else
    {
      fi[3] = sqp0p3;
      fi[4] = fi[1];
    }
    fi[5] = fi[1];
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL
  void imzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fi = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fi[0] = cxmake( pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf );
    fi[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv chi = cxmake( -(fptype)nhel * fpsqrt( -2. * pvec3 ), 0. );
    fi[3] = cxzero_sv();
    fi[4] = cxzero_sv();
    if ( nh == 1 )
    {
      fi[2] = cxzero_sv();
      fi[5] = chi;
    }
    else
    {
      fi[2] = chi;
      fi[5] = cxzero_sv();
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL
  void ixzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fi = reinterpret_cast<cxtype_sv*>( wavefunctions );
    //fi[0] = cxmake( -pvec0 * nsf, -pvec2 * nsf ); // AV: BUG! not the same as ixxxxx
    //fi[1] = cxmake( -pvec0 * nsf, -pvec1 * nsf ); // AV: BUG! not the same as ixxxxx
    fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf ); // AV: BUG FIX
    fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf ); // AV: BUG FIX
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1/sqp0p3, pvec2/sqp0p3 );
    if ( nh == 1 )
    {
      fi[2] = cxzero_sv();
      fi[3] = cxzero_sv();
      fi[4] = chi0;
      fi[5] = chi1;
    }
    else
    {
      fi[2] = chi1;
      fi[3] = chi0;
      fi[4] = cxzero_sv();
      fi[5] = cxzero_sv();
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void vxxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* vc = reinterpret_cast<cxtype_sv*>( wavefunctions );
    const fptype sqh = fpsqrt( 0.5 ); // AV this is > 0!
    const fptype hel = nhel;
    vc[0] = cxmake( pvec0 * (fptype)nsv, pvec3 * (fptype)nsv );
    vc[1] = cxmake( pvec1 * (fptype)nsv, pvec2 * (fptype)nsv );
    if ( vmass != 0. )
    {
      const int nsvahl = nsv * sycl::fabs( hel );
      const fptype_sv pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
      const fptype_sv pp = fpmin( pvec0, fpsqrt( pt2 + ( pvec3 * pvec3 ) ) );
      const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
      const fptype hel0 = 1. - sycl::fabs( hel );
      if ( pp == 0. )
      {
        vc[2] = cxmake( 0., 0. );
        vc[3] = cxmake( -hel * sqh, 0. );
        vc[4] = cxmake( 0., nsvahl * sqh );
        vc[5] = cxmake( hel0, 0. );
      }
      else
      {
        const fptype emp = pvec0 / ( vmass * pp );
        vc[2] = cxmake( hel0 * pp / vmass, 0. );
        vc[5] = cxmake( hel0 * pvec3 * emp + hel * pt / pp * sqh, 0. );
        if ( pt != 0. )
        {
          const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
          vc[3] = cxmake( hel0 * pvec1 * emp - pvec1 * pzpt, -nsvahl * pvec2 / pt * sqh );
          vc[4] = cxmake( hel0 * pvec2 * emp - pvec2 * pzpt, nsvahl * pvec1 / pt * sqh );
        }
        else
        {
          vc[3] = cxmake( -hel * sqh, 0. );
          vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
        }
      }
    }
    else
    {
      const fptype_sv& pp = pvec0; // NB: rewrite the following as in Fortran, using pp instead of pvec0
      const fptype_sv pt = fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) );
      vc[2] = cxzero_sv();
      vc[5] = cxmake( hel * pt / pp * sqh, 0. );
      if ( pt != 0. )
      {
        const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
        vc[3] = cxmake( -pvec1 * pzpt, -nsv * pvec2 / pt * sqh );
        vc[4] = cxmake( -pvec2 * pzpt, nsv * pvec1 / pt * sqh );
      }
      else
      {
        vc[3] = cxmake( -hel * sqh, 0. );
        vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void sxxxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype,                 // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               //const int,                    // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* sc = reinterpret_cast<cxtype_sv*>( wavefunctions );
    sc[2] = cxmake( 1 + fptype_sv{0}, 0 );
    sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
    sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void oxxxxx( const fptype momenta[],
               const size_t ievt,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fo = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    if ( fmass != 0. )
    {
      const fptype_sv pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
      if ( pp == 0. )
      {
        fptype sqm[2] = { fpsqrt( sycl::fabs( fmass ) ), 0. }; // possibility of negative fermion masses
        //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
        sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
        const int ip = -( ( 1 - nh ) / 2 ) * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        const int im = ( 1 + nh ) / 2 * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
        fo[2] = cxmake( im * sqm[sycl::abs( ip )], 0 );
        fo[3] = cxmake( ip * nsf * sqm[sycl::abs( ip )], 0 );
        fo[4] = cxmake( im * nsf * sqm[sycl::abs( im )], 0 );
        fo[5] = cxmake( ip * sqm[sycl::abs( im )], 0 );
      }
      else
      {
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype omega[2] = { fpsqrt( pvec0 + pp ), 0. };
        omega[1] = fmass / omega[0];
        const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
        const fptype sfomeg[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype pp3 = fpmax( pp + pvec3, 0. );
        const cxtype chi[2] = { cxmake( fpsqrt( pp3 * (fptype)0.5 / pp ), 0. ),
                                ( ( pp3 == 0. ) ? cxmake( -nh, 0. )
                                  : cxmake( nh * pvec1, -pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        fo[2] = sfomeg[1] * chi[im];
        fo[3] = sfomeg[1] * chi[ip];
        fo[4] = sfomeg[0] * chi[im];
        fo[5] = sfomeg[0] * chi[ip];
      }
    }
    else
    {
      const fptype_sv sqp0p3 = fpternary( ( pvec1 == 0. ) and ( pvec2 == 0. ) and ( pvec3 < 0. ),
                                          0, fpsqrt( fpmax( pvec0 + pvec3, 0. ) ) * (fptype)nsf );
      const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ),
                                 cxternary( ( sqp0p3 == 0. ),
                                            cxmake( -nhel, 0. ) * fpsqrt( 2. * pvec0 ),
                                            cxmake( (fptype)nh * pvec1, -pvec2 ) / sqp0p3 ) };
      if ( nh == 1 )
      {
        fo[2] = chi[0];
        fo[3] = chi[1];
        fo[4] = cxzero_sv();
        fo[5] = cxzero_sv();
      }
      else
      {
        fo[2] = cxzero_sv();
        fo[3] = cxzero_sv();
        fo[4] = chi[1];
        fo[5] = chi[0];
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL
  void opzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fo = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fo[0] = cxmake( pvec3 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv csqp0p3 = cxmake( fpsqrt( 2. * pvec3 ) * (fptype)nsf, 0. );
    fo[3] = cxzero_sv();
    fo[4] = cxzero_sv();
    if ( nh == 1 )
    {
      fo[2] = csqp0p3;
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[5] = csqp0p3;
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL
  void omzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fo = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fo[0] = cxmake( -pvec3 * (fptype)nsf, pvec3 * (fptype)nsf ); // remember pvec0 == -pvec3
    fo[1] = cxzero_sv();
    const int nh = nhel * nsf;
    const cxtype_sv chi1 = cxmake( -nhel, 0. ) * fpsqrt( -2. * pvec3 );
    if ( nh == 1 )
    {
      fo[2] = cxzero_sv();
      fo[3] = chi1;
      fo[4] = cxzero_sv();
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[3] = cxzero_sv();
      fo[4] = chi1;
      //fo[5] = chi1; // AV: BUG!
      fo[5] = cxzero_sv(); // AV: BUG FIX
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL
  void oxzxxx( const fptype momenta[],
               const size_t ievt,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               fptype wavefunctions[],
               const int ipar )                // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& pvec0 = kernelAccessIp4IparIevt( momenta, 0, ipar, ievt );
    const fptype_sv& pvec1 = kernelAccessIp4IparIevt( momenta, 1, ipar, ievt );
    const fptype_sv& pvec2 = kernelAccessIp4IparIevt( momenta, 2, ipar, ievt );
    const fptype_sv& pvec3 = kernelAccessIp4IparIevt( momenta, 3, ipar, ievt );
    cxtype_sv* fo = reinterpret_cast<cxtype_sv*>( wavefunctions );
    fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
    fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
    const int nh = nhel * nsf;
    //const float sqp0p3 = sqrtf( pvec0 + pvec3 ) * nsf; // AV: why force a float here?
    const fptype_sv sqp0p3 = fpsqrt( pvec0 + pvec3 ) * (fptype)nsf;
    const cxtype_sv chi0 = cxmake( sqp0p3, 0. );
    const cxtype_sv chi1 = cxmake( (fptype)nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
    if ( nh == 1 )
    {
      fo[2] = chi0;
      fo[3] = chi1;
      fo[4] = cxzero_sv();
      fo[5] = cxzero_sv();
    }
    else
    {
      fo[2] = cxzero_sv();
      fo[3] = cxzero_sv();
      fo[4] = chi1;
      fo[5] = chi0;
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //==========================================================================
