

template<typename T>
KOKKOS_FUNCTION void ixxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex_t<double> fi[6])
{
  
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + fmass * fmass);
  fi[0] = complex_t<double> (-p0 * nsf, -pvec(3) * nsf);
  fi[1] = complex_t<double> (-pvec(1) * nsf, -pvec(2) * nsf);
  if (fmass != 0.0)
  {
    double pp = min(p0, sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3)));
    if (pp == 0.0)
    {
      double sqm[2] = {sqrt(std::abs(fmass)), 0.};
      sqm[1] = (fmass < 0) ? -sqm[0] : sqm[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      fi[2] = complex_t<double>(ip * sqm[ip], 0.);
      fi[3] = complex_t<double>(im * nsf * sqm[ip], 0.);
      fi[4] = complex_t<double>(ip * nsf * sqm[im], 0.);
      fi[5] = complex_t<double>(im * sqm[im], 0.);
    }
    else
    {
      const double sf[2] = { (1 + nsf + (1 - nsf) * (nhel * nsf)) * 0.5,
                             (1 + nsf - (1 - nsf) * (nhel * nsf)) * 0.5 };
      double omega[2] = { sqrt(p0 + pp), 0. };
      omega[1] = fmass / omega[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      const double sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const double pp3 = max(pp + pvec(3), 0.0);
      
      const complex_t<double> chi[2] = { complex_t<double> (sqrt(pp3 * 0.5/pp), 0),
                                               ( pp3 == 0. ?
                                                 complex_t<double> (-(nhel * nsf), 0) :
                                                 complex_t<double> ((nhel * nsf) * pvec(1), pvec(2)) 
                                                                  / sqrt(2.0 * pp * pp3) ) };
      fi[2] = sfomega[0] * chi[im];
      fi[3] = sfomega[0] * chi[ip];
      fi[4] = sfomega[1] * chi[im];
      fi[5] = sfomega[1] * chi[ip];
    }
  }
  else
  { 
    const double sqp0p3 = (pvec(1) == 0.0 and pvec(2) == 0.0 and pvec(3) < 0.0) ? 0. : sqrt(max(p0 + pvec(3), 0.0)) * nsf;

    const complex_t<double> chi[2] = {complex_t<double> (sqp0p3, 0.0),
          (sqp0p3 == 0.0) ? complex_t<double> (-nhel * sqrt(2.0 * p0), 0.0) :
                            complex_t<double> ((nhel * nsf) * pvec(1), pvec(2))/sqp0p3 };
    if ((nhel * nsf) == 1)
    {
      fi[2] = complex_t<double> (0.0, 0.0);
      fi[3] = complex_t<double> (0.0, 0.0);
      fi[4] = chi[0];
      fi[5] = chi[1];
    }
    else
    {
      fi[2] = chi[1];
      fi[3] = chi[0];
      fi[4] = complex_t<double> (0.0, 0.0);
      fi[5] = complex_t<double> (0.0, 0.0);
    }
  }
  return;
}


template<typename T>
KOKKOS_FUNCTION void ipzxxx(const T& pvec, const int& nhel, const int& nsf, complex_t<double> fi[6])
{
  const double& pvec3 = pvec(3);
  
  fi[0] = complex_t<double> (-pvec3 * nsf, -pvec3 * nsf);
  fi[1] = complex_t<double> (0.,0.);

  complex_t<double> sqp0p3 = complex_t<double>(sqrt(2.* pvec3) * nsf, 0.);
  const int nh = nhel * nsf;
  fi[2]=fi[1];
  if(nh==1){
   fi[3] = fi[1];
   fi[4] = sqp0p3;
  }else{
   fi[3] = sqp0p3;
   fi[4] = fi[1];
  }
  fi[5]=fi[1];
}


template<typename T>
KOKKOS_FUNCTION void imzxxx(const T& pvec, const int nhel, const int nsf, complex_t<double> fi[6])
{

  const double& pvec3 = pvec(3);
  fi[0] = complex_t<double> (pvec3 * nsf, -pvec3 * nsf);
  fi[1] = complex_t<double> (0., 0.);
  complex_t<double>  chi = complex_t<double> (-nhel * sqrt(-2.0 * pvec3), 0.0);
  const int nh = nhel * nsf;
  fi[3]=fi[1];
  fi[4]=fi[1];
  if (nh ==1) {
     fi[2] = fi[1];
     fi[5] = chi;
  }else{
     fi[2] = chi;
     fi[5] = fi[1];
  }    
}


template<typename T>
KOKKOS_FUNCTION void ixzxxx(const T& pvec, const int& nhel, const int& nsf, complex_t<double> fi[6])
{
  const double& pvec0 = pvec(0);
  const double& pvec1 = pvec(1);
  const double& pvec2 = pvec(2);
  const double& pvec3 = pvec(3);

  fi[0] = complex_t<double> (-pvec0 * nsf, -pvec2 * nsf);
  fi[1] = complex_t<double> (-pvec0 * nsf, -pvec1 * nsf);
  const int nh = nhel * nsf;

  float sqp0p3 = sqrtf(pvec0 + pvec3) * nsf;
  complex_t<double> chi0 = complex_t<double> (sqp0p3, 0.0);
  complex_t<double> chi1 = complex_t<double> (nh * pvec1/sqp0p3, pvec2/sqp0p3);
  complex_t<double> CZERO = complex_t<double>(0.,0.);

  if (nh ==1){
    fi[2]=CZERO;
    fi[3]=CZERO;
    fi[4]=chi0 ;
    fi[5]=chi1 ;
  }else{
     fi[2]= chi1;
     fi[3]= chi0;
     fi[4]= CZERO;
     fi[5]= CZERO;
  }

  return;
}


template<typename T>
KOKKOS_FUNCTION void vxxxxx(const T& pvec, const double vmass, 
  const int nhel, const int nsv, complex_t<double> vc[6]) 
{
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + vmass * vmass);
  const double sqh = sqrt(0.5);
  const double hel = nhel;
  const int nsvahl = nsv * std::abs(hel);
  vc[0] = complex_t<double> ( p0 * nsv, pvec(3) * nsv );
  vc[1] = complex_t<double> ( pvec(1) * nsv, pvec(2) * nsv );
  if ( vmass != 0. )
  {
    const double pt2 = (pvec(1) * pvec(1)) + (pvec(2) * pvec(2));
    const double pp = std::min(p0, sqrt(pt2 + (pvec(3) * pvec(3))));
    const double pt = std::min(pp, sqrt(pt2));
    const double hel0 = 1. - std::abs( hel );
    if ( pp == 0. )
    {
      vc[2] = complex_t<double> ( 0., 0. );
      vc[3] = complex_t<double> ( -hel * sqh, 0. );
      vc[4] = complex_t<double> ( 0., nsvahl * sqh );
      vc[5] = complex_t<double> ( hel0, 0. );
    }
    else
    {
      const double emp = p0 / ( vmass * pp );
      vc[2] = complex_t<double> ( hel0 * pp / vmass, 0. );
      vc[5] = 
      complex_t<double> ( hel0 * pvec(3) * emp + hel * pt / pp * sqh, 0. );
      if ( pt != 0. )
      {
        const double pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
        vc[3] = complex_t<double> ( hel0 * pvec(1) * emp - pvec(1) * pzpt, -nsvahl * pvec(2) / pt * sqh);
        vc[4] = complex_t<double> ( hel0 * pvec(2) * emp - pvec(2) * pzpt, nsvahl * pvec(1) / pt * sqh);
      }
      else
      {
        vc[3] = complex_t<double> ( -hel * sqh, 0. );
        vc[4] = complex_t<double> (0., nsvahl * (pvec(3) < 0  ? -sqh : sqh ) );
      }
    }
  }
  else
  {
    const double pt = sqrt((pvec(1) * pvec(1)) + (pvec(2) * pvec(2)));
    vc[2] = complex_t<double> (0., 0. );
    vc[5] = complex_t<double> (hel * pt / p0 * sqh, 0. );
    if ( pt != 0. )
    {
      const double pzpt = pvec(3) / (p0 * pt) * sqh * hel;
      vc[3] = complex_t<double> ( -pvec(1) * pzpt, -nsv * pvec(2) / pt * sqh);
      vc[4] = complex_t<double> ( -pvec(2) * pzpt, nsv * pvec(1) / pt * sqh);
    }
    else
    {
      vc[3] = complex_t<double> ( -hel * sqh, 0. );
      vc[4] = complex_t<double> ( 0., nsv * ( pvec(3) < 0 ? -sqh : sqh ) );
    }
  }
  return;
}


template<typename T>
KOKKOS_FUNCTION  void sxxxxx(const T& pvec, const double& smass, const int& nhel, const int& nss, complex_t<double> sc[3])
{
  const double& p0 = pvec(0);
  const double& p1 = pvec(1);
  const double& p2 = pvec(2);
  const double& p3 = pvec(3);
  //double p[4] = {0, pvec[0], pvec[1], pvec[2]};
  //p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+smass*smass);
  sc[2] = complex_t<double>(1.00, 0.00);
  sc[0] = complex_t<double>(p0 * nss, p3 * nss);
  sc[1] = complex_t<double>(p1 * nss, p2 * nss);
  return;
}


template<typename T>
KOKKOS_FUNCTION void oxxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex_t<double> fo[6]) 
{
  
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + fmass * fmass);

  fo[0] = complex_t<double> (p0 * nsf, pvec(3) * nsf);
  fo[1] = complex_t<double> (pvec(1) * nsf, pvec(2) * nsf);
  if (fmass != 0.)
  {
    const double pp = min(p0, sqrt((pvec(1) * pvec(1)) + (pvec(2) * pvec(2)) + (pvec(3) * pvec(3))));
    if (pp == 0.)
    {
      double sqm[2] = {sqrt(std::abs(fmass)), 0.};
      sqm[1] = fmass < 0 ? -sqm[0] : sqm[0];
      const int ip = -((1 - (nhel * nsf))/2) * nhel;
      const int im = (1 + (nhel * nsf))/2 * nhel;
      fo[2] = complex_t<double> (im * sqm[std::abs(ip)], 0.);
      fo[3] = complex_t<double> (ip * nsf * sqm[std::abs(ip)], 0.);
      fo[4] = complex_t<double> (im * nsf * sqm[std::abs(im)], 0.);
      fo[5] = complex_t<double> (ip * sqm[std::abs(im)], 0.);
    }
    else
    {
      const double sf[]= {double(1 + nsf + (1 - nsf) * (nhel * nsf)) * 0.5,
                          double(1 + nsf - (1 - nsf) * (nhel * nsf)) * 0.5};
      double omega[2] = {sqrt(p0 + pp), 0.};
      omega[1] = fmass/omega[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      double sfomeg[2] = { sf[0] * omega[ip], sf[1] * omega[im]};
      const double pp3 = max(pp + pvec(3), 0.);
      const complex_t<double> chi[2] = { complex_t<double> (sqrt(pp3 * 0.5/pp), 0.),
                    pp3 == 0. ? complex_t<double> (-(nhel * nsf), 0.00) :
                                complex_t<double> ((nhel * nsf) * pvec(1), -pvec(2))/sqrt(2.0 * pp * pp3)
        };
      fo[2] = sfomeg[1] * chi[im];
      fo[3] = sfomeg[1] * chi[ip];
      fo[4] = sfomeg[0] * chi[im];
      fo[5] = sfomeg[0] * chi[ip];
    }
  }
  else
  {
    double sqp0p3 = ((pvec(1) == 0.00) and (pvec(2) == 0.00) and (pvec(3) < 0.00)) ?
              0. : sqrt(max(p0 + pvec(3), 0.00)) * nsf;

    const complex_t<double> chi[2] = { complex_t<double> (sqp0p3, 0.00),
            sqp0p3 == 0. ? complex_t<double> (-nhel, 0.) * sqrt(2. * p0) :
                           complex_t<double> ((nhel * nsf) * pvec(1), -pvec(2))/sqp0p3
    };

    if ((nhel * nsf) == 1)
    {
      fo[2] = chi[0];
      fo[3] = chi[1];
      fo[4] = complex_t<double> (0., 0.);
      fo[5] = complex_t<double> (0., 0.);
    }
    else
    {
      fo[2] = complex_t<double> (0., 0.);
      fo[3] = complex_t<double> (0., 0.);
      fo[4] = chi[1];
      fo[5] = chi[0];
    }
  }
  return;
}

template<typename T>
KOKKOS_FUNCTION  void opzxxx(const T& pvec, const int& nhel, const int& nsf, complex_t<double> fo[6])
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = PZ

  const double& pvec3 = pvec(3);

  fo[0] = complex_t<double> (pvec3 * nsf, pvec3 * nsf);
  fo[1] = complex_t<double> (0., 0.);
  const int nh = nhel * nsf;

  complex_t<double> CSQP0P3 = complex_t<double> (sqrt(2.* pvec3) * nsf, 0.00);


  fo[3]=fo[1];
  fo[4]=fo[1];
  if (nh==1){
     fo[2]=CSQP0P3;
     fo[5]=fo[1];
  }else{
     fo[2]=fo[1];
     fo[5]=CSQP0P3;
  }
}


template<typename T>
KOKKOS_FUNCTION  void omzxxx(const T& pvec, const int& nhel, const int& nsf, complex_t<double> fo[6])
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = -PZ (E>0)

  const double& pvec3 = pvec(3);
  fo[0] = complex_t<double> (-pvec3 * nsf, pvec3 * nsf);
  fo[1] = complex_t<double> (0., 0.);
  const int nh = nhel * nsf;
  complex_t<double> chi = complex_t<double> (-nhel, 0.00) * sqrt(-2.0 * pvec3);

  if(nh==1){
    fo[2]=fo[1];
    fo[3]=chi;
    fo[4]=fo[1];
    fo[5]=fo[1];
  }else{
    fo[2]= fo[1];
    fo[3]= fo[1];
    fo[4]= chi;
    fo[5]= chi;
  }
  return;
}

template<typename T>
KOKKOS_FUNCTION  void oxzxxx(const T& pvec, const int& nhel, const int& nsf, complex_t<double> fo[6])
{
  // ASSUMPTIONS FMASS =0
  // PT > 0
  const double& p0 = pvec(0);
  const double& p1 = pvec(1);
  const double& p2 = pvec(2);
  const double& p3 = pvec(3);

  fo[0] = complex_t<double> (p0 * nsf, p3 * nsf);
  fo[1] = complex_t<double> (p1 * nsf, p2 * nsf);
  const int nh = nhel * nsf;

  float sqp0p3 = sqrtf(p0 + p3) * nsf;
  complex_t<double> chi0 = complex_t<double> (sqp0p3, 0.00);
  complex_t<double> chi1 = complex_t<double> (nh * p1/sqp0p3, -p2/sqp0p3);
  complex_t<double> zero = complex_t<double> (0.00, 0.00);

  if(nh==1){
    fo[2]=chi0;
    fo[3]=chi1;
    fo[4]=zero;
    fo[5]=zero;
  }else{
    fo[2]=zero;
    fo[3]=zero;
    fo[4]=chi1;
    fo[5]=chi0;
  }
  return;
}

