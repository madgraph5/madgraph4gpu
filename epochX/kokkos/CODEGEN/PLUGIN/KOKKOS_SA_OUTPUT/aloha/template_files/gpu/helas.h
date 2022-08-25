
  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void ixxxxx(const T& pvec,
                            const fptype fmass, 
                            const int nhel,
                            const int nsf,
                            cxtype_sv* fi)
{
  
  
  fi[0] = cxmake(-pvec(0) * nsf, -pvec(3) * nsf);
  fi[1] = cxmake(-pvec(1) * nsf, -pvec(2) * nsf);
  const int nh = nhel * nsf;
  if ( fmass != 0. )
  {
    const fptype_sv pp = fpmin(pvec(0), fpsqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) ) );
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
      fptype omega[2] = { fpsqrt( pvec(0) + pp ), 0. };
      omega[1] = fmass / omega[0];
      const int ip = ( 1 + nh ) / 2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
      const int im = ( 1 - nh ) / 2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
      const fptype sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const fptype pp3 = fpmax( pp + pvec(3), 0. );
      const cxtype chi[2] = { cxmake( fpsqrt ( pp3 * (fptype)0.5 / pp ), 0. ),
                              ( pp3 == 0. ?
                                cxmake( -nh, 0. ) :
                                cxmake( nh * pvec(1), pvec(2) ) / fpsqrt( 2. * pp * pp3 ) ) };
      fi[2] = sfomega[0] * chi[im];
      fi[3] = sfomega[0] * chi[ip];
      fi[4] = sfomega[1] * chi[im];
      fi[5] = sfomega[1] * chi[ip];
    }
    
  }
  else
  { 
    const fptype_sv sqp0p3 = fpternary( ( pvec(1) == 0. and pvec(2) == 0. and pvec(3) < 0. ),
                                        fptype_sv{0}, fpsqrt(fpmax(pvec(0) + pvec(3), 0. ) ) * (fptype)nsf );

    const cxtype_sv chi[2] = { cxmake( sqp0p3, 0. ), cxternary( ( sqp0p3 == 0. ),
                                                               cxmake( -(fptype)nhel * fpsqrt( 2. * pvec(0) ), 0. ),
                                                               cxmake( (fptype)nh * pvec(1), pvec(2)) / sqp0p3 ) };
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void ipzxxx(const T& pvec,
                            const int& nhel,
                            const int& nsf,
                            cxtype_sv fi[])
{
  const fptype& pvec3 = pvec(3);
  
  fi[0] = cxmake(-pvec3 * (fptype)nsf, -pvec3 * (fptype)nsf);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void imzxxx(const T& pvec, 
                            const int nhel,
                            const int nsf,
                            cxtype_sv fi[])
{

  const fptype pvec3 = pvec(3);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION void ixzxxx(const T& pvec,
                            const int& nhel,
                            const int& nsf,
                            cxtype_sv fi[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);

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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void vxxxxx(const T& pvec,
                            const fptype vmass,
                            const int nhel,
                            const int nsv,
                            cxtype_sv vc[]) 
{
  const fptype sqh = fpsqrt( 0.5 );
  const fptype hel = nhel;
  const int nsvahl = nsv * abs(hel);
  vc[0] = cxmake( pvec(0) * (fptype)nsv, pvec(3) * (fptype)nsv );
  vc[1] = cxmake( pvec(1) * (fptype)nsv, pvec(2) * (fptype)nsv );
  if ( vmass != 0. )
  {
    const int nsvahl = nsv * std::abs( hel );
    const fptype_sv pt2 = ( pvec(1) * pvec(1) ) + ( pvec(2) * pvec(2) );
    const fptype_sv pp = fpmin( pvec(0), fpsqrt( pt2 + ( pvec(3) * pvec(3) ) ) );
    const fptype_sv pt = fpmin( pp, fpsqrt( pt2 ) );
    const fptype hel0 = 1. - std::abs( hel );
    if ( pp == 0. )
    {
      vc[2] = cxmake( 0., 0. );
      vc[3] = cxmake( -hel * sqh, 0. );
      vc[4] = cxmake( 0., nsvahl * sqh );
      vc[5] = cxmake( hel0, 0. );
    }
    else
    {
      const fptype emp = pvec(0) / ( vmass * pp );
      vc[2] = cxmake( hel0 * pp / vmass, 0. );
      vc[5] = cxmake( hel0 * pvec(3) * emp + hel * pt / pp * sqh, 0. );
      if ( pt != 0. )
      {
        const fptype pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
        vc[3] = cxmake( hel0 * pvec(1) * emp - pvec(1) * pzpt, -nsvahl * pvec(2) / pt * sqh );
        vc[4] = cxmake( hel0 * pvec(2) * emp - pvec(2) * pzpt, nsvahl * pvec(1) / pt * sqh );
      }
      else
      {
        vc[3] = cxmake( -hel * sqh, 0. );
            // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
            //vc[4] = cxmake( 0., nsvahl * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV: why abs here?
        vc[4] = cxmake( 0., nsvahl * (pvec(3) < 0 ? -sqh : sqh ) );
      }
    }
  }
  else
  {
    const fptype_sv& pp = pvec(0); // NB: rewrite the following as in Fortran, using pp instead of pvec0
    const fptype_sv pt = fpsqrt( ( pvec(1) * pvec(1) ) + ( pvec(2) * pvec(2) ) );
    vc[2] = cxzero_sv();
    vc[5] = cxmake( hel * pt / pp * sqh, 0. );
    if ( pt != 0. )
    {
      const fptype pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
      vc[3] = cxmake( -pvec(1) * pzpt, -nsv * pvec(2) / pt * sqh );
      vc[4] = cxmake( -pvec(2) * pzpt, nsv * pvec(1) / pt * sqh );
    }
    else
    {
      vc[3] = cxmake( -hel * sqh, 0. );
          // NB: Do not use "abs" for floats! It returns an integer with no build warning! Use std::abs!
          //vc[4] = cxmake( 0, nsv * ( pvec3 < 0. ? -std::abs( sqh ) : std::abs( sqh ) ) ); // AV why abs here?
      vc[4] = cxmake( 0., nsv * ( pvec(3) < 0. ? -sqh : sqh ) ); // AV: removed an abs here
    }
  }
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION  void sxxxxx(const T& pvec,
                             const fptype& smass,
                             const int& nhel,
                             const int& nss,
                             cxtype_sv sc[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  sc[2] = cxmake( 1 + fptype_sv{0}, 0 );
      sc[0] = cxmake( pvec0 * (fptype)nss, pvec3 * (fptype)nss );
      sc[1] = cxmake( pvec1 * (fptype)nss, pvec2 * (fptype)nss );
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
template<typename T>
KOKKOS_INLINE_FUNCTION void oxxxxx(const T& pvec,
                            const fptype fmass,
                            const int nhel,
                            const int nsf,
                            cxtype_sv fo[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
  fo[0] = cxmake( pvec0 * (fptype)nsf, pvec3 * (fptype)nsf );
  fo[1] = cxmake( pvec1 * (fptype)nsf, pvec2 * (fptype)nsf );
  const int nh = nhel * nsf;
  if ( fmass != 0. )
  {
    const fptype_sv pp = fpmin( pvec0, fpsqrt( ( pvec1 * pvec1 ) + ( pvec2 * pvec2 ) + ( pvec3 * pvec3 ) ) );
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void opzxxx(const T& pvec,
                             const int& nhel,                  // input: -1 or +1 (helicity of fermion)
                             const int& nsf,                  // input: +1 (particle) or -1 (antiparticle)
                             cxtype_sv fo[])
{

  const fptype pvec3 = pvec(3);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void omzxxx(const T& pvec,
                             const int& nhel,
                             const int& nsf,
                             cxtype_sv fo[])
{

  const fptype pvec3 = pvec(3);
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
  return;
}

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
template<typename T>
KOKKOS_INLINE_FUNCTION  void oxzxxx(const T& pvec,
                             const int& nhel,
                             const int& nsf,
                             cxtype_sv fo[])
{
  const fptype pvec0 = pvec(0);
  const fptype pvec1 = pvec(1);
  const fptype pvec2 = pvec(2);
  const fptype pvec3 = pvec(3);
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
  return;
}


