  __device__
  inline const fptype& pIparIp4Ievt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    const int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    //printf( "%f\n", momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

#ifndef __CUDACC__
  // Return by value: it seems a tiny bit faster than returning a reference (both for scalar and vector), not clear why
  inline fptype_sv pIparIp4Ipag( const fptype_sv* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                 const int ipar,
                                 const int ip4,
                                 const int ipagM )
  {
#ifndef MGONGPU_CPPSIMD
    // NB THERE IS NO SIMD YET IN GGTTGG! HENCE ipagM=ievt
    return pIparIp4Ievt( momenta1d, ipar, ip4, ipagM );
#else
#error THERE IS NO SIMD YET IN GGTTGG
#endif
  }
#endif
  
  //--------------------------------------------------------------------------

  __device__
  void ixxxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly slower, in c++
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      //const fptype pvec0 = fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ); // AV: BUG?! (NOT AS IN THE FORTRAN)
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt ); // AV: BUG FIX (DO AS IN THE FORTRAN)
#else
      //printf( "ixxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      fi[0] = cxmake( -pvec0 * (fptype)nsf, -pvec3 * (fptype)nsf );
      fi[1] = cxmake( -pvec1 * (fptype)nsf, -pvec2 * (fptype)nsf );
      const int nh = nhel * nsf;
      if ( fmass != 0. )
      {
        const fptype_sv pp = fpmin( pvec0, fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ) );
#ifndef MGONGPU_CPPSIMD
        if ( pp == 0. )
        {
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
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
#else
        const int ip = ( 1 + nh ) / 2;
        const int im = ( 1 - nh ) / 2;
        // Branch A: pp == 0.
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0 }; // possibility of negative fermion masses (NB: SCALAR!)
        sqm[1] = ( fmass < 0 ? -sqm[0] : sqm[0] ); // AV: removed an abs here (as above)
        const cxtype fiA_2 = ip * sqm[ip]; // scalar cxtype: real part initialised from fptype, imag part = 0
        const cxtype fiA_3 = im * nsf * sqm[ip]; // scalar cxtype: real part initialised from fptype, imag part = 0
        const cxtype fiA_4 = ip * nsf * sqm[im]; // scalar cxtype: real part initialised from fptype, imag part = 0
        const cxtype fiA_5 = im * sqm[im]; // scalar cxtype: real part initialised from fptype, imag part = 0
        // Branch B: pp != 0.
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_v omega[2] = { fpsqrt( pvec0 + pp ), 0 };
        omega[1] = fmass / omega[0];
        const fptype_v sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
        const fptype_v pp3 = fpmax( pp + pvec3, 0 );
        const cxtype_v chi[2] = { cxmake( fpsqrt ( pp3 * 0.5 / pp ), 0 ),
                                  cxternary( ( pp3 == 0. ),
                                             cxmake( -nh, 0 ),
                                             cxmake( (fptype)nh * pvec1, pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) };
        const cxtype_v fiB_2 = sfomega[0] * chi[im];
        const cxtype_v fiB_3 = sfomega[0] * chi[ip];
        const cxtype_v fiB_4 = sfomega[1] * chi[im];
        const cxtype_v fiB_5 = sfomega[1] * chi[ip];
        // Choose between the results from branch A and branch B
        const bool_v mask = ( pp == 0. );
        fi[2] = cxternary( mask, fiA_2, fiB_2 );
        fi[3] = cxternary( mask, fiA_3, fiB_3 );
        fi[4] = cxternary( mask, fiA_4, fiB_4 );
        fi[5] = cxternary( mask, fiA_5, fiB_5 );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void ipzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ipzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "ipzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void imzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "imzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void ixzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "ixzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void vxxxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,              // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,               // input: +1 (final) or -1 (initial)
               cxtype_sv* vc,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "vxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "vxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      const fptype sqh = fpsqrt( 0.5 ); // AV this is > 0!
      const fptype hel = nhel;
      vc[0] = cxmake( pvec0 * (fptype)nsv, pvec3 * (fptype)nsv );
      vc[1] = cxmake( pvec1 * (fptype)nsv, pvec2 * (fptype)nsv );
      if ( vmass != 0. )
      {
        const int nsvahl = nsv * std::abs( hel );
        const fptype_sv pt2 = ( pvec1 * pvec1 ) + ( pvec2 * pvec2 );
        const fptype_sv pp = fpmin( pvec0, fpsqrt( pt2 + ( pvec3 * pvec3 ) ) );
        const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
        const fptype hel0 = 1. - std::abs( hel );
#ifndef MGONGPU_CPPSIMD
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
            // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
            //vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
            vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
          }
        }
#else
        // Branch A: pp == 0.
        const cxtype vcA_2 = cxmake( 0, 0 );
        const cxtype vcA_3 = cxmake( -hel * sqh, 0 );
        const cxtype vcA_4 = cxmake( 0, nsvahl * sqh );
        const cxtype vcA_5 = cxmake( hel0, 0 );
        // Branch B: pp != 0.
        const fptype_v emp = pvec0 / ( vmass * pp );
        const cxtype_v vcB_2 = cxmake( hel0 * pp / vmass, 0 );
        const cxtype_v vcB_5 = cxmake( hel0 * pvec3 * emp + hel * pt / pp * sqh, 0 );
        // Branch B1: pp != 0. and pt != 0.
        const fptype_v pzpt = pvec3 / ( pp * pt ) * sqh * hel;
        const cxtype_v vcB1_3 = cxmake( hel0 * pvec1 * emp - pvec1 * pzpt, -(fptype)nsvahl * pvec2 / pt * sqh );
        const cxtype_v vcB1_4 = cxmake( hel0 * pvec2 * emp - pvec2 * pzpt, (fptype)nsvahl * pvec1 / pt * sqh );
        // Branch B2: pp != 0. and pt == 0.
        const cxtype vcB2_3 = cxmake( -hel * sqh, 0. );
        const cxtype_v vcB2_4 = cxmake( 0., (fptype)nsvahl * fpternary( ( pvec3 < 0 ), -sqh, sqh ) ); // AV: removed an abs here
        // Choose between the results from branch A and branch B (and from branch B1 and branch B2)
        const bool_v mask = ( pp == 0. );
        const bool_v maskB = ( pt != 0. );
        vc[2] = cxternary( mask, vcA_2, vcB_2 );
        vc[3] = cxternary( mask, vcA_3, cxternary( maskB, vcB1_3, vcB2_3 ) );
        vc[4] = cxternary( mask, vcA_4, cxternary( maskB, vcB1_4, vcB2_4 ) );
        vc[5] = cxternary( mask, vcA_5, vcB_5 );
#endif
      }
      else
      {
        const fptype_sv& pp = pvec0; // NB: rewrite the following  as in Fortran, using pp instead of pvec0
        const fptype_sv pt = fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) );
        vc[2] = cxzero_sv();
        vc[5] = cxmake( hel * pt / pp * sqh, 0. );
#ifndef MGONGPU_CPPSIMD
        if ( pt != 0. )
        {
          const fptype pzpt = pvec3 / ( pp * pt ) * sqh * hel;
          vc[3] = cxmake( -pvec1 * pzpt, -nsv * pvec2 / pt * sqh );
          vc[4] = cxmake( -pvec2 * pzpt, nsv * pvec1 / pt * sqh );
        }
        else
        {
          vc[3] = cxmake( -hel * sqh, 0. );
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          //vc[4] = cxmake( 0, nsv * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV why abs here?
          vc[4] = cxmake( 0., nsv * ( pvec3 < 0. ? -sqh : sqh ) ); // AV: removed an abs here
        }
#else
        // Branch A: pt != 0.
        const fptype_v pzpt = pvec3 / ( pp * pt ) * sqh * hel;
        const cxtype_v vcA_3 = cxmake( -pvec1 * pzpt, -(fptype)nsv * pvec2 / pt * sqh );
        const cxtype_v vcA_4 = cxmake( -pvec2 * pzpt, (fptype)nsv * pvec1 / pt * sqh );
        // Branch B: pt == 0.
        const cxtype vcB_3 = cxmake( -(fptype)hel * sqh, 0 );
        const cxtype_v vcB_4 = cxmake( 0, (fptype)nsv * fpternary( ( pvec3 < 0 ), -sqh, sqh ) ); // AV: removed an abs here
        // Choose between the results from branch A and branch B
        const bool_v mask = ( pt != 0. );
        vc[3] = cxternary( mask, vcA_3, vcB_3 );
        vc[4] = cxternary( mask, vcA_4, vcB_4 );
#endif
      }
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void sxxxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype,                // WARNING: "smass" unused (missing in Fortran)
               const int,                   // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,               // input: +1 (final) or -1 (initial)
               cxtype_sv sc[3],             // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "sxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "sxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      sc[2] = cxmake( 1 + fptype_sv{0}, 0 );
      sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
      sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void oxxxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly faster, in c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "oxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
      fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
      const int nh = nhel * nsf;
      if ( fmass != 0. )
      {
        const fptype_sv pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
#ifndef MGONGPU_CPPSIMD
        if ( pp == 0. )
        {
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0. }; // possibility of negative fermion masses
          //sqm[1] = ( fmass < 0. ? -abs( sqm[0] ) : abs( sqm[0] ) ); // AV: why abs here?
          sqm[1] = ( fmass < 0. ? -sqm[0] : sqm[0] ); // AV: removed an abs here
          const int ip = -( ( 1 - nh ) / 2 ) * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          const int im = ( 1 + nh ) / 2 * nhel; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          fo[2] = cxmake( im * sqm[std::abs( ip )], 0 );
          fo[3] = cxmake( ip * nsf * sqm[std::abs( ip )], 0 );
          fo[4] = cxmake( im * nsf * sqm[std::abs( im )], 0 );
          fo[5] = cxmake( ip * sqm[std::abs( im )], 0 );
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
#else
        // Branch A: pp == 0.
        // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
        fptype sqm[2] = { fpsqrt( std::abs( fmass ) ), 0 }; // possibility of negative fermion masses
        sqm[1] = ( fmass < 0 ? -sqm[0] : sqm[0] ); // AV: removed an abs here (as above)
        const int ipA = -( ( 1 - nh ) / 2 ) * nhel;
        const int imA = ( 1 + nh ) / 2 * nhel;
        const cxtype foA_2 = imA * sqm[std::abs( ipA )];
        const cxtype foA_3 = ipA * nsf * sqm[std::abs( ipA )];
        const cxtype foA_4 = imA * nsf * sqm[std::abs( imA )];
        const cxtype foA_5 = ipA * sqm[std::abs( imA )];
        // Branch B: pp != 0.
        const fptype sf[2] = { fptype( 1 + nsf + ( 1 - nsf ) * nh ) * (fptype)0.5,
                               fptype( 1 + nsf - ( 1 - nsf ) * nh ) * (fptype)0.5 };
        fptype_v omega[2] = { fpsqrt( pvec0 + pp ), 0 };
        omega[1] = fmass / omega[0];
        const int ipB = ( 1 + nh ) / 2;
        const int imB = ( 1 - nh ) / 2;
        const fptype_v sfomeg[2] = { sf[0] * omega[ipB], sf[1] * omega[imB] };
        const fptype_v pp3 = fpmax( pp + pvec3, 0. );
        const cxtype_v chi[2] = { cxmake( fpsqrt( pp3 * 0.5 / pp ), 0. ),
                                  ( cxternary( ( pp3 == 0. ),
                                               cxmake( -nh, 0. ),
                                               cxmake( (fptype)nh * pvec1, -pvec2 ) / fpsqrt( 2. * pp * pp3 ) ) ) };
        const cxtype_v foB_2 = sfomeg[1] * chi[imB];
        const cxtype_v foB_3 = sfomeg[1] * chi[ipB];
        const cxtype_v foB_4 = sfomeg[0] * chi[imB];
        const cxtype_v foB_5 = sfomeg[0] * chi[ipB];
        // Choose between the results from branch A and branch B
        const bool_v mask = ( pp == 0. );
        fo[2] = cxternary( mask, foA_2, foB_2 );
        fo[3] = cxternary( mask, foA_3, foB_3 );
        fo[4] = cxternary( mask, foA_4, foB_4 );
        fo[5] = cxternary( mask, foA_5, foB_5 );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void opzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "opzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "opzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void omzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ipzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "ipzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void oxzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const int nhel,              // input: -1 or +1 (helicity of fermion)
               const int nsf,               // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "oxzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
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
    }
    // +++ END EVENT LOOP (where necessary) +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

//--------------------------------------------------------------------------
/*
__device__
inline const fptype& pIparIp4Ievt(const fptype * momenta1d,  // input: momenta as AOSOA[npagM][npar][4][neppM]
const int ipar, 
const int ip4, 
const int ievt)
{
  // mapping for the various scheme AOS, OSA, ...

  using mgOnGpu::np4; 
  using mgOnGpu::npar; 
  const int neppM = mgOnGpu::neppM;  // ASA layout: constant at compile-time
  fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM])
      momenta1d; // cast to multiD array pointer (AOSOA)
  const int ipagM = ievt/neppM;  // #eventpage in this iteration
  const int ieppM = ievt%neppM;  // #event in the current eventpage in this iteration
  // return allmomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + ip4*neppM +
  // ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  return momenta[ipagM][ipar][ip4][ieppM]; 
}

//--------------------------------------------------------------------------

__device__ void ixxxxx(const fptype * allmomenta, const fptype& fmass, const
    int& nhel, const int& nsf,
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  mgDebug(0, __FUNCTION__); 
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
#ifndef __CUDACC__
  using std::max; 
  using std::min; 
#endif

  const fptype& pvec0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
  const fptype& pvec1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& pvec2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  cxtype chi[2]; 
  fptype sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 

  fptype p[4] = {0, pvec1, pvec2, pvec3}; 
  //p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass); // AV: BUG?! (NOT AS IN THE FORTRAN)
  p[0] = pvec0; // AV: BUG FIX (DO AS IN THE FORTRAN)
  fi[0] = cxtype(-p[0] * nsf, -p[3] * nsf); 
  fi[1] = cxtype(-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])); 
    if (pp == 0.0)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      //sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]); // BUG issue #198
      sqm[1] = (fmass < 0) ? -sqm[0] : sqm[0];
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = cxtype(sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = cxtype(-nh, 0); 
      }
      else
      {
        chi[1] = 
        cxtype(nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf; 
    }
    chi[0] = cxtype(sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = cxtype(-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = cxtype(nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = cxtype(0.0, 0.0); 
      fi[3] = cxtype(0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = cxtype(0.0, 0.0); 
      fi[5] = cxtype(0.0, 0.0); 
    }
  }
  //++ END LOOP ON IEVT ++
  mgDebug(1, __FUNCTION__); 
  return; 
}


__device__ void ipzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTION FMASS == 0
  // PX = PY = 0
  // E = P3 (E>0)
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  // const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  fi[0] = cxtype (-pvec3 * nsf, -pvec3 * nsf); 
  fi[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 

  cxtype sqp0p3 = cxtype(sqrt(2. * pvec3) * nsf, 0.); 

  fi[2] = fi[1]; 
  if(nh == 1)
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
}

__device__ void imzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTION FMASS == 0
  // PX = PY = 0
  // E = -P3 (E>0)
  // printf("p3 %f", pvec[2]);
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  fi[0] = cxtype (pvec3 * nsf, -pvec3 * nsf); 
  fi[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 
  cxtype chi = cxtype (-nhel * sqrt(-2.0 * pvec3), 0.0); 


  fi[3] = fi[1]; 
  fi[4] = fi[1]; 
  if (nh == 1)
  {
    fi[2] = fi[1]; 
    fi[5] = chi; 
  }
  else
  {
    fi[2] = chi; 
    fi[5] = fi[1]; 
  }
}

__device__ void ixzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS: FMASS == 0
  // Px and Py are not zero

  // cxtype chi[2];
  // fptype sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2];
  // int ip, im, nh;
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& pvec1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& pvec2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 


  // fptype p[4] = {(float), (float) pvec[0], (float) pvec[1], (float) pvec[2]};
  // p[0] = sqrtf(p[3] * p[3] + p[1] * p[1] + p[2] * p[2]);

  fi[0] = cxtype (-pvec0 * nsf, -pvec2 * nsf); 
  fi[1] = cxtype (-pvec0 * nsf, -pvec1 * nsf); 
  int nh = nhel * nsf; 

  float sqp0p3 = sqrtf(pvec0 + pvec3) * nsf; 
  cxtype chi0 = cxtype (sqp0p3, 0.0); 
  cxtype chi1 = cxtype (nh * pvec1/sqp0p3, pvec2/sqp0p3); 
  cxtype CZERO = cxtype(0., 0.); 

  if (nh == 1)
  {
    fi[2] = CZERO; 
    fi[3] = CZERO; 
    fi[4] = chi0; 
    fi[5] = chi1; 
  }
  else
  {
    fi[2] = chi1; 
    fi[3] = chi0; 
    fi[4] = CZERO; 
    fi[5] = CZERO; 
  }

  return; 
}

__device__ void vxxxxx(const fptype * allmomenta, const fptype& vmass, const
    int& nhel, const int& nsv,
cxtype vc[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  fptype hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 

#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#else
  using std::min; 
#endif

  const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+vmass*vmass);

  sqh = sqrt(0.5); 
  hel = fptype(nhel); 
  nsvahl = nsv * std::abs(hel); 
  pt2 = (p1 * p1) + (p2 * p2); 
  pp = min(p0, sqrt(pt2 + (p3 * p3))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = cxtype(p0 * nsv, p3 * nsv); 
  vc[1] = cxtype(p1 * nsv, p2 * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - std::abs(hel); 
    if (pp == 0.0)
    {
      vc[2] = cxtype(0.0, 0.0); 
      vc[3] = cxtype(-hel * sqh, 0.0); 
      vc[4] = cxtype(0.0, nsvahl * sqh); 
      vc[5] = cxtype(hel0, 0.0); 
    }
    else
    {
      emp = p0/(vmass * pp); 
      vc[2] = cxtype(hel0 * pp/vmass, 0.0); 
      vc[5] = 
      cxtype(hel0 * p3 * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p3/(pp * pt) * sqh * hel; 
        vc[3] = cxtype(hel0 * p1 * emp - p1 * pzpt, 
         - nsvahl * p2/pt * sqh); 
        vc[4] = cxtype(hel0 * p2 * emp - p2 * pzpt, 
        nsvahl * p1/pt * sqh); 
      }
      else
      {
        vc[3] = cxtype(-hel * sqh, 0.0); 
        //vc[4] = cxtype(0.0, nsvahl * (p3 < 0) ? - abs(sqh) : abs(sqh)); // BUG issue #198
        vc[4] = cxtype(0.0, nsvahl * (p3 < 0) ? -sqh : sqh);
      }
    }
  }
  else
  {
    // pp = p0;
    pt = sqrt((p1 * p1) + (p2 * p2)); 
    vc[2] = cxtype(0.0, 0.0); 
    vc[5] = cxtype(hel * pt/p0 * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p3/(p0 * pt) * sqh * hel; 
      vc[3] = cxtype(-p1 * pzpt, -nsv * p2/pt * sqh); 
      vc[4] = cxtype(-p2 * pzpt, nsv * p1/pt * sqh); 
    }
    else
    {
      vc[3] = cxtype(-hel * sqh, 0.0); 
      //vc[4] = cxtype(0.0, nsv * (p3 < 0) ? - abs(sqh) : abs(sqh)); // BUG issue #198
      vc[4] = cxtype(0.0, nsv * (p3 < 0) ? -sqh : sqh);
    }
  }
  return; 
}

__device__ void sxxxxx(const fptype * allmomenta, const fptype& smass, const
    int& nhel, const int& nss,
cxtype sc[3], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)
{
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
  const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+smass*smass);
  sc[2] = cxtype(1.00, 0.00); 
  sc[0] = cxtype(p0 * nss, p3 * nss); 
  sc[1] = cxtype(p1 * nss, p2 * nss); 
  return; 
}

__device__ void oxxxxx(const fptype * allmomenta, const fptype& fmass, const
    int& nhel, const int& nsf,
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
#ifndef __CUDACC__
  using std::min; 
  using std::max; 
#endif
  cxtype chi[2]; 
  fptype sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 

  const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+fmass*fmass);

  fo[0] = cxtype(p0 * nsf, p3 * nsf); 
  fo[1] = cxtype(p1 * nsf, p2 * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p0, sqrt((p1 * p1) + (p2 * p2) + (p3 * p3))); 
    if (pp == 0.000)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      //sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]); // BUG issue #198
      sqm[1] = (fmass < 0) ? -sqm[0] : sqm[0];
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[std::abs(ip)]; 
      fo[3] = ip * nsf * sqm[std::abs(ip)]; 
      fo[4] = im * nsf * sqm[std::abs(im)]; 
      fo[5] = ip * sqm[std::abs(im)]; 
    }
    else
    {
      sf[0] = fptype(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = fptype(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p0 + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p3, 0.00); 
      chi[0] = cxtype(sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = cxtype(-nh, 0.00); 
      }
      else
      {
        chi[1] = 
        cxtype(nh * p1, -p2)/sqrt(2.0 * pp * pp3); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if ((p1 == 0.00) and (p2 == 0.00) and (p3 < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = sqrt(max(p0 + p3, 0.00)) * nsf; 
    }
    chi[0] = cxtype(sqp0p3, 0.00); 
    if (sqp0p3 == 0.000)
    {
      chi[1] = cxtype(-nhel, 0.00) * sqrt(2.0 * p0); 
    }
    else
    {
      chi[1] = cxtype(nh * p1, -p2)/sqp0p3; 
    }
    if (nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = cxtype(0.00, 0.00); 
      fo[5] = cxtype(0.00, 0.00); 
    }
    else
    {
      fo[2] = cxtype(0.00, 0.00); 
      fo[3] = cxtype(0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

__device__ void opzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = PZ
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  fo[0] = cxtype (pvec3 * nsf, pvec3 * nsf); 
  fo[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 

  cxtype CSQP0P3 = cxtype (sqrt(2. * pvec3) * nsf, 0.00); 


  fo[3] = fo[1]; 
  fo[4] = fo[1]; 
  if (nh == 1)
  {
    fo[2] = CSQP0P3; 
    fo[5] = fo[1]; 
  }
  else
  {
    fo[2] = fo[1]; 
    fo[5] = CSQP0P3; 
  }
}


__device__ void omzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = -PZ (E>0)
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  fo[0] = cxtype (-pvec3 * nsf, pvec3 * nsf); 
  fo[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 
  cxtype chi = cxtype (-nhel, 0.00) * sqrt(-2.0 * pvec3); 

  if(nh == 1)
  {
    fo[2] = fo[1]; 
    fo[3] = chi; 
    fo[4] = fo[1]; 
    fo[5] = fo[1]; 
  }
  else
  {
    fo[2] = fo[1]; 
    fo[3] = fo[1]; 
    fo[4] = chi; 
    fo[5] = chi; 
  }
  return; 
}

__device__ void oxzxxx(const fptype * allmomenta, const int& nhel, const int&
    nsf,
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PT > 0
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  // float p[4] = {0, (float) pvec[0], (float) pvec[1], (float) pvec[2]};
  // p[0] = sqrtf(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);

  fo[0] = cxtype (p0 * nsf, p3 * nsf); 
  fo[1] = cxtype (p1 * nsf, p2 * nsf); 
  int nh = nhel * nsf; 

  float sqp0p3 = sqrtf(p0 + p3) * nsf; 
  cxtype chi0 = cxtype (sqp0p3, 0.00); 
  cxtype chi1 = cxtype (nh * p1/sqp0p3, -p2/sqp0p3); 
  cxtype zero = cxtype (0.00, 0.00); 

  if(nh == 1)
  {
    fo[2] = chi0; 
    fo[3] = chi1; 
    fo[4] = zero; 
    fo[5] = zero; 
  }
  else
  {
    fo[2] = zero; 
    fo[3] = zero; 
    fo[4] = chi1; 
    fo[5] = chi0; 
  }
  return; 
}
*/
//--------------------------------------------------------------------------

__device__ constexpr fptype one( 1. );
