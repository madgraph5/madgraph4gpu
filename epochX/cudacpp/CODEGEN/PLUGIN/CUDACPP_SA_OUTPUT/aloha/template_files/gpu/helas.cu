  //--------------------------------------------------------------------------

  // Decode momentum AOSOA: compute address of fptype for the given particle, 4-momentum component and event
  // Return the fptype by reference (equivalent to returning its memory address)
  __device__ inline
  const fptype& pIparIp4Ievt( const fptype* momenta, // input: momenta as AOSOA[npagM][npar][4][neppM]
                              const int ipar,
                              const int ip4,
                              const int ievt )
  {
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    //printf( "%2d %2d %8d %8.3f\n", ipar, ip4, ievt, momenta[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return momenta[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
    //fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta; // cast to multiD array pointer (AOSOA)
    //return momenta[ipagM][ipar][ip4][ieppM]; // this seems ~1-2% faster in eemumu C++?
  }

#ifndef __CUDACC__
  // Return a SIMD vector of fptype's for neppV events (for the given particle, 4-momentum component and event "V-page")
  // Return the vector by value: strictly speaking, this is only unavoidable for neppM<neppV
  // For neppM>=neppV (both being powers of 2), the momentum neppM-AOSOA is reinterpreted in terms of neppV-vectors:
  // it could also be returned by reference, but no performance degradation is observed when returning by value
  inline
  fptype_sv pIparIp4Ipag( const fptype* momenta, // input: momenta as AOSOA[npagM][npar][4][neppM]
                          const int ipar,
                          const int ip4,
                          const int ipagV )
  {
    const int ievt0 = ipagV*neppV; // virtual event V-page ipagV contains neppV events [ievt0...ievt0+neppV-1]
#ifdef MGONGPU_CPPSIMD
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    // Use c++17 "if constexpr": compile-time branching
    if constexpr ( ( neppM >= neppV ) && ( neppM%neppV == 0 ) )
    {
      constexpr bool useReinterpretCastIfPossible = true; // DEFAULT
      //constexpr bool useReinterpretCastIfPossible = false; // FOR PERFORMANCE TESTS
      constexpr bool skipAlignmentCheck = true; // DEFAULT (MAY SEGFAULT, NEEDS A SANITY CHECK ELSEWHERE!)
      //constexpr bool skipAlignmentCheck = false; // SLOWER BUT SAFER
      if constexpr ( useReinterpretCastIfPossible && skipAlignmentCheck )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! skip alignment check" << std::endl; first=false; } // SLOWS DOWN...
        // Fastest (4.92E6 in eemumu 512y)
        // This requires alignment for momenta - segmentation fault otherwise!
        return *reinterpret_cast<const fptype_sv*>( &( pIparIp4Ievt( momenta, ipar, ip4, ievt0 ) ) );
      }
      else if ( useReinterpretCastIfPossible && ( (size_t)(momenta) % mgOnGpu::cppAlign == 0 ) )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! alignment ok, use reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (6%) slower (4.62E6 in eemumu 512y) because of the alignment check
        // This explicitly checks alignment for momenta to avoid segmentation faults
        return *reinterpret_cast<const fptype_sv*>( &( pIparIp4Ievt( momenta, ipar, ip4, ievt0 ) ) );
      }
      else
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! AOSOA but no reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (2%) slower (4.86E6 in eemumu 512y)
        // This does not require alignment for momenta, but it requires AOSOA with neppM>=neppV and neppM%neppV==0
        const fptype* out0 = &( pIparIp4Ievt( momenta, ipar, ip4, ievt0) );
#if MGONGPU_CPPSIMD == 2
        return fptype_v{ *( out0 ),
                         *( out0+1 ) };
#elif MGONGPU_CPPSIMD == 4
        return fptype_v{ *( out0 ),
                         *( out0+1 ),
                         *( out0+2 ),
                         *( out0+3 ) };
#elif MGONGPU_CPPSIMD == 8
        return fptype_v{ *( out0 ),
                         *( out0+1 ),
                         *( out0+2 ),
                         *( out0+3 ),
                         *( out0+4 ),
                         *( out0+5 ),
                         *( out0+6 ),
                         *( out0+7 ) };
#elif MGONGPU_CPPSIMD == 16
        return fptype_v{ *( out0 ),
                         *( out0+1 ),
                         *( out0+2 ),
                         *( out0+3 ),
                         *( out0+4 ),
                         *( out0+5 ),
                         *( out0+6 ),
                         *( out0+7 ),
                         *( out0+8 ),
                         *( out0+9 ),
                         *( out0+10 ),
                         *( out0+11 ),
                         *( out0+12 ),
                         *( out0+13 ),
                         *( out0+14 ),
                         *( out0+15 ) };
#else
#warning Internal error? This code should never be reached!
        // A bit (2%) slower (4.86E6 in eemumu 512y)
        // This does not require alignment for momenta, but it requires AOSOA with neppM>=neppV and neppM%neppV==0
        fptype_v out;
        for ( int ieppV=0; ieppV<neppV; ieppV++ ) out[ieppV] = *( out0 + ieppV );
        return out;
#endif
      }
    }
    else
    {
      // Much (20%) slower (4.07E6 in eemumu 512y)
      // This does not even require AOSOA with neppM>=neppV and neppM%neppV==0 (e.g. can be used with AOS neppM==1)
#if MGONGPU_CPPSIMD == 2
      return fptype_v{ pIparIp4Ievt( momenta, ipar, ip4, ievt0 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+1 ) };
#elif MGONGPU_CPPSIMD == 4
      return fptype_v{ pIparIp4Ievt( momenta, ipar, ip4, ievt0 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+1 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+2 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+3 ) };
#elif MGONGPU_CPPSIMD == 8
      return fptype_v{ pIparIp4Ievt( momenta, ipar, ip4, ievt0 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+1 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+2 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+3 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+4 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+5 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+6 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+7 ) };
#elif MGONGPU_CPPSIMD == 16
      return fptype_v{ pIparIp4Ievt( momenta, ipar, ip4, ievt0 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+1 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+2 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+3 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+4 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+5 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+6 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+7 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+8 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+9 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+10 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+11 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+12 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+13 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+14 ),
                       pIparIp4Ievt( momenta, ipar, ip4, ievt0+15 ) };
#else
#warning Internal error? This code should never be reached!
      // Much much (40%) slower (3.02E6 in eemumu 512y)
      fptype_v out;
      for ( int ieppV=0; ieppV<neppV; ieppV++ ) out[ieppV] = pIparIp4Ievt( momenta, ipar, ip4, ievt0+ieppV );
      return out;
#endif
    }
#else
    return pIparIp4Ievt( momenta, ipar, ip4, ievt0 );
#endif
  }
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  __device__
  void ixxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly slower, in c++
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
      //const fptype pvec0 = fpsqrt( pvec1 * pvec1 + pvec2 * pvec2 + pvec3 * pvec3 ); // AV: BUG?! (NOT AS IN THE FORTRAN)
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt ); // AV: BUG FIX (DO AS IN THE FORTRAN)
#else
      //printf( "ixxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__
  void ipzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ipzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "ipzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__
  void imzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "imzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__
  void ixzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "ixzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  __device__
  void vxxxxx( const fptype* momenta,
               const fptype vmass,             // input: vector boson mass
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv vc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "vxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "vxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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
        const fptype_sv& pp = pvec0; // NB: rewrite the following as in Fortran, using pp instead of pvec0
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

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  __device__
  void sxxxxx( const fptype* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const int,                      // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "sxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "sxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  __device__
  void oxxxxx( const fptype* momenta,
               const fptype fmass,             // input: fermion mass
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxxxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems faster in cuda, in spite of more registers used
      // AV: copying by value (not by ref) seems irrelevant, or slightly faster, in c++
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "oxxxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  __device__
  void opzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "opzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "opzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  __device__
  void omzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ipzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copy by value (not by ref) as this seems faster in cuda for other functions
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "ipzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  __device__
  void oxzxxx( const fptype* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[],
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
    // +++ START EVENT LOOP (where necessary) +++
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      // AV: copying by value (not by ref) seems to give the same performance in both cuda and c++
      const fptype pvec0 = pIparIp4Ievt( momenta, ipar, 0, ievt );
      const fptype pvec1 = pIparIp4Ievt( momenta, ipar, 1, ievt );
      const fptype pvec2 = pIparIp4Ievt( momenta, ipar, 2, ievt );
      const fptype pvec3 = pIparIp4Ievt( momenta, ipar, 3, ievt );
#else
      //printf( "oxzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( momenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( momenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( momenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( momenta, ipar, 3, ipagV );
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
