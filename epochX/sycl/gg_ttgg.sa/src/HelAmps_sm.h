// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Sep 2010) for the MG5aMC backend.
//==========================================================================
// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
//==========================================================================
// Copyright (C) 2021-2023 Argonne National Laboratory.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
//==========================================================================
// This file has been automatically generated for SYCL standalone by
// MadGraph5_aMC@NLO v. 3.5.0_lo_vect, 2023-01-26
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

//#include <cmath>
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm
{
! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Sep 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2023 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: O. Mattelaer (Mar 2020) for the MG5aMC CUDACPP plugin.
! Further modified by: O. Mattelaer, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
!==========================================================================
! Copyright (C) 2021-2023 Argonne National Laboratory.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
!==========================================================================
  //--------------------------------------------------------------------------

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__( ( always_inline ) )
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void ixxxxx( const vector4* momenta,
               const fptype fmass,             // input: fermion mass
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL INLINE
  void ipzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL INLINE
  void imzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL INLINE
  void ixzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void vxxxxx( const vector4* momenta,
               const fptype vmass,             // input: vector boson mass
               const signed char nhel,         // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const signed char nsv,          // input: +1 (final) or -1 (initial)
               cxtype_sv vc[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void sxxxxx( const vector4* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const signed char nhel,         // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const signed char nss,          // input: +1 (final) or -1 (initial)
               cxtype_sv sc[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void oxxxxx( const vector4* momenta,
               const fptype fmass,             // input: fermion mass
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL INLINE
  void opzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL INLINE
  void omzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL INLINE
  void oxzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void ixxxxx( const vector4* momenta,
               const fptype fmass,             // input: fermion mass
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;

      fi[0] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec0, fptype(-nsf)*pvec3);
      fi[1] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec1, fptype(-nsf)*pvec2);
      const signed char nh = nhel*nsf;
      if (fmass != 0.0) {
          auto pp = FPMIN(pvec0, FPSQRT(pvec1*pvec1 + pvec2*pvec2 + pvec3*pvec3));
          auto pp_zero = pp == 0.0; 

          fptype sqm[2];
          sqm[0] = FPSQRT(sycl::fabs(fmass)); // possibility of negative fermion masses
          sqm[1] = fmass < 0.0 ? -sqm[0] : sqm[0];

          const signed char ip = (1 + nh)/2;
          const signed char im = (1 - nh)/2;

          const fptype sf[2] = { 0.5*(1 + nsf + (1 - nsf)*nh),
                                 0.5*(1 + nsf - (1 - nsf)*nh)  };
          fptype_sv omega[2];
          omega[0] = FPSQRT(pvec0 + pp);
          omega[1] = fmass/omega[0];
          const fptype_sv sfomega[2] = { sf[0]*omega[ip], sf[1]*omega[im] };
          auto pp3 = FPMAX(pp + pvec3, 0.0);
          auto pp3_zero = pp3 == 0.0;
          const cxtype_sv chi[2] = { CXMAKE_SV_2ARG(FPSQRT(0.5*pp3/pp), FPZERO_SV),
                                     CXMAKE_SV_2ARG(FPCONDITIONAL_SV(fptype(nh)*pvec1/FPSQRT(2.0*pp*pp3), fptype_sv(-nh), pp3_zero),
                                                    FPCONDITIONAL_SV(pvec2/FPSQRT(2.0*pp*pp3), FPZERO_SV, pp3_zero))
                                   };
          fi[2] = CXCONDITIONAL_SV(sfomega[0]*chi[im], CXMAKE_SV_2ARG(fptype_sv(ip*sqm[ip]), FPZERO_SV)     , pp_zero);
          fi[3] = CXCONDITIONAL_SV(sfomega[0]*chi[ip], CXMAKE_SV_2ARG(fptype_sv(im*nsf*sqm[ip]), FPZERO_SV) , pp_zero);
          fi[4] = CXCONDITIONAL_SV(sfomega[1]*chi[im], CXMAKE_SV_2ARG(fptype_sv(ip*nsf*sqm[im]), FPZERO_SV) , pp_zero);
          fi[5] = CXCONDITIONAL_SV(sfomega[1]*chi[ip], CXMAKE_SV_2ARG(fptype_sv(im*sqm[im]), FPZERO_SV)     , pp_zero);
      }
      else {
          const fptype_sv sqp0p3 = FPCONDITIONAL_SV(fptype(nsf)*FPSQRT(FPMAX(pvec0 + pvec3, FPZERO_SV)), pvec1, (pvec1 == 0.0) && (pvec2 == 0.0) && (pvec3 < 0.0));
          const cxtype_sv chi[2] = { CXMAKE_SV_2ARG(sqp0p3, FPZERO_SV),
                                     CXCONDITIONAL_SV(CXMAKE_SV_2ARG(fptype(nh)*pvec1/sqp0p3, pvec2/sqp0p3),
                                                      CXMAKE_SV_2ARG(fptype(-nhel)*FPSQRT(2.0*pvec0), FPZERO_SV),
                                                      sqp0p3 == 0.0) };
          if (nh == 1) {
              fi[2] = CXZERO_SV;
              fi[3] = CXZERO_SV;
              fi[4] = chi[0];
              fi[5] = chi[1];
          }
          else {
              fi[2] = chi[1];
              fi[3] = chi[0];
              fi[4] = CXZERO_SV;
              fi[5] = CXZERO_SV;
          }
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL
  void ipzxxx( const vector4* momenta,
               //const fptype fmass,                   // ASSUME fermion mass==0
               const signed char nhel,                 // input: -1 or +1 (helicity of fermion)
               const signed char nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) {
      auto pvec3 = momenta[0].z;
      fi[0] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec3, fptype(-nsf)*pvec3);
      fi[1] = CXZERO_SV;
      const cxtype_sv sqp0p3 = CXMAKE_SV_2ARG(fptype(nsf)*FPSQRT(2.0*pvec3), FPZERO_SV);
      fi[2] = fi[1];
      if (nhel*nsf == 1) {
          fi[3] = fi[1];
          fi[4] = sqp0p3;
      }
      else {
          fi[3] = sqp0p3;
          fi[4] = fi[1];
      }
      fi[5] = fi[1];
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL
  void imzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) {
      auto pvec3 = momenta[0].z;
      fi[0] = CXMAKE_SV_2ARG(fptype(nsf)*pvec3, fptype(-nsf)*pvec3);
      fi[1] = CXZERO_SV;
      const cxtype_sv chi = CXMAKE_SV_2ARG(fptype(-nhel)*FPSQRT(-2.0*pvec3), FPZERO_SV);
      fi[3] = CXZERO_SV;
      fi[4] = CXZERO_SV;
      if (nhel*nsf == 1) {
          fi[2] = CXZERO_SV;
          fi[5] = chi;
      }
      else {
          fi[2] = chi;
          fi[5] = CXZERO_SV;
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fi[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL
  void ixzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fi[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;
      fi[0] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec0, fptype(-nsf)*pvec3);
      fi[1] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec1, fptype(-nsf)*pvec2);
      const fptype_sv sqp0p3 = fptype(nsf)*FPSQRT(pvec0 + pvec3);
      const cxtype_sv chi0 = CXMAKE_SV_2ARG(sqp0p3, FPZERO_SV);
      const cxtype_sv chi1 = CXMAKE_SV_2ARG(fptype(nhel*nsf)*pvec1/sqp0p3, pvec2/sqp0p3);
      if (nhel*nsf == 1) {
          fi[2] = CXZERO_SV;
          fi[3] = CXZERO_SV;
          fi[4] = chi0;
          fi[5] = chi1;
      }
      else {
          fi[2] = chi1;
          fi[3] = chi0;
          fi[4] = CXZERO_SV;
          fi[5] = CXZERO_SV;
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction vc[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void vxxxxx( const vector4* momenta,
               const fptype vmass,             // input: vector boson mass
               const signed char nhel,         // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const signed char nsv,          // input: +1 (final) or -1 (initial)
               cxtype_sv vc[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;
      vc[0] = CXMAKE_SV_2ARG(fptype(nsv)*pvec0, fptype(nsv)*pvec3);
      vc[1] = CXMAKE_SV_2ARG(fptype(nsv)*pvec1, fptype(nsv)*pvec2);
      if (vmass != 0.0) {
          const signed char nsvahl = nsv*sycl::abs(nhel);
          const fptype_sv pt2 = pvec1*pvec1 + pvec2*pvec2;
          const fptype_sv pp = FPMIN(pvec0, FPSQRT(pt2 + pvec3*pvec3));
          const fptype_sv pt = FPMIN(pp, FPSQRT(pt2));
          const signed char hel0 = 1 - sycl::abs(nhel);
          auto pp_zero = pp == 0.0; 
          auto pt_zero = pt == 0.0; 
          auto emp = pvec0/pp/vmass;
          auto pzpt = fptype(nhel)*SQRTH*pvec3/pp/pt;

          vc[2] = CXMAKE_SV_2ARG(FPCONDITIONAL_SV(hel0/vmass*pp, FPZERO_SV, pp_zero), FPZERO_SV);
          vc[3] = CXCONDITIONAL_SV(CXCONDITIONAL_SV(CXMAKE_SV_2ARG(fptype(hel0)*pvec1*emp - pvec1*pzpt, fptype(-nsvahl)*SQRTH*pvec2/pt),
                                                    CXMAKE_SV_2ARG(fptype_sv(fptype(-nhel)*SQRTH), FPZERO_SV), pt_zero),
                               CXMAKE_SV_2ARG(fptype_sv(fptype(-nhel)*SQRTH), FPZERO_SV), pp_zero);
          vc[4] = CXCONDITIONAL_SV(CXCONDITIONAL_SV(CXMAKE_SV_2ARG(fptype(hel0)*pvec2*emp - pvec2*pzpt, fptype(nsvahl)*SQRTH*pvec1/pt),
                                                    CXMAKE_SV_2ARG(FPZERO_SV, FPCONDITIONAL_SV(fptype_sv(fptype(nsvahl)*SQRTH), fptype_sv(fptype(-nsvahl)*SQRTH), pvec3 < 0.0)), pt_zero),
                               CXMAKE_SV_2ARG(FPZERO_SV, fptype_sv(fptype(nsvahl)*SQRTH)), pp_zero);
          vc[5] = CXMAKE_SV_2ARG(FPCONDITIONAL_SV(fptype(hel0)*pvec3*emp + fptype(nhel)*SQRTH*pt/pp, fptype_sv(fptype(hel0)), pp_zero), FPZERO_SV);
      }
      else {
          const fptype_sv& pp = pvec0; // NB: rewrite the following as in Fortran, using pp instead of pvec0
          const fptype_sv pt = FPSQRT(pvec1*pvec1 + pvec2*pvec2);
          auto pt_zero = pt == 0.0;
          auto pzpt = fptype(nhel)*SQRTH*pvec3/pp/pt;

          vc[2] = CXZERO_SV;
          vc[3] = CXCONDITIONAL_SV(CXMAKE_SV_2ARG(-pvec1*pzpt, fptype(-nsv)*SQRTH*pvec2/pt), CXMAKE_SV_2ARG(fptype_sv(fptype(-nhel)*SQRTH), FPZERO_SV), pt_zero);
          vc[4] = CXCONDITIONAL_SV(CXMAKE_SV_2ARG(-pvec2*pzpt, fptype(nsv)*SQRTH*pvec1/pt),
                                   CXMAKE_SV_2ARG(FPZERO_SV, FPCONDITIONAL_SV(fptype_sv(nsv*SQRTH), fptype_sv(fptype(-nsv)*SQRTH), pvec3 < 0.0)), pt_zero);
          vc[5] = CXMAKE_SV_2ARG(fptype(nhel)*SQRTH*pt/pp, FPZERO_SV);
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction sc[3] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void sxxxxx( const vector4* momenta,
               const fptype,                   // WARNING: input "smass" unused (missing in Fortran) - scalar boson mass
               const signed char nhel,         // WARNING: input "nhel" unused (missing in Fortran) - scalar has no helicity!
               const signed char nss,          // input: +1 (final) or -1 (initial)
               cxtype_sv sc[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;

      sc[2] = CXMAKE_SV_2ARG(FPONE_SV, FPZERO_SV);
      sc[0] = CXMAKE_SV_2ARG(fptype(nss)*pvec0, fptype(nss)*pvec3);
      sc[1] = CXMAKE_SV_2ARG(fptype(nss)*pvec1, fptype(nss)*pvec2);
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  SYCL_EXTERNAL
  void oxxxxx( const vector4* momenta,
               const fptype fmass,             // input: fermion mass
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;

      fo[0] = CXMAKE_SV_2ARG(fptype(nsf)*pvec0, fptype(nsf)*pvec3);
      fo[1] = CXMAKE_SV_2ARG(fptype(nsf)*pvec1, fptype(nsf)*pvec2);
      const signed char nh = nhel*nsf;
      if (fmass != 0.0) {
          const fptype_sv pp = FPMIN(pvec0, FPSQRT(pvec1*pvec1 + pvec2*pvec2 + pvec3*pvec3));
          auto pp_zero = pp == 0.0; 
          fptype sqm[2];
          sqm[0] = FPSQRT(sycl::fabs(fmass)); // possibility of negative fermion masses
          sqm[1] = fmass < 0.0 ? -sqm[0] : sqm[0];
          const signed char ip_zero = nhel*(nh - 1)/2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++
          const signed char im_zero = nhel*(1 + nh)/2; // NB: Fortran sqm(0:1) also has indexes 0,1 as in C++

          const fptype sf[2] = { 0.5*(1 + nsf + (1 - nsf)*nh),
                                 0.5*(1 + nsf - (1 - nsf)*nh)  };
          fptype_sv omega[2];
          omega[0] = FPSQRT(pvec0 + pp);
          omega[1] = fmass/omega[0];
          const signed char ip = (1 + nh)/2; // NB: Fortran is (3+nh)/2 because omega(2) has indexes 1,2 and not 0,1
          const signed char im = (1 - nh)/2; // NB: Fortran is (3-nh)/2 because omega(2) has indexes 1,2 and not 0,1
          const fptype_sv sfomeg[2] = { sf[0]*omega[ip], sf[1]*omega[im] };
          const fptype_sv pp3 = FPMAX(pp + pvec3, 0.0);
          const cxtype_sv chi[2] = { CXMAKE_SV_2ARG(FPSQRT(0.5*pp3/pp), FPZERO_SV),
                                  CXCONDITIONAL_SV(CXMAKE_SV_2ARG(fptype(nh)*pvec1/FPSQRT(2.0*pp*pp3), -pvec2/FPSQRT(2.0*pp*pp3)), CXMAKE_SV_2ARG(fptype_sv(fptype(-nh)), FPZERO_SV), pp3 == 0.0)
                                }; 

          fo[2] = CXCONDITIONAL_SV(sfomeg[1]*chi[im], CXMAKE_SV_2ARG(fptype_sv(fptype(im_zero)*sqm[sycl::abs(ip_zero)]), FPZERO_SV), pp_zero);
          fo[3] = CXCONDITIONAL_SV(sfomeg[1]*chi[ip], CXMAKE_SV_2ARG(fptype_sv(fptype(ip_zero)*nsf*sqm[sycl::abs(ip_zero)]), FPZERO_SV), pp_zero);
          fo[4] = CXCONDITIONAL_SV(sfomeg[0]*chi[im], CXMAKE_SV_2ARG(fptype_sv(fptype(im_zero)*nsf*sqm[sycl::abs(im_zero)]), FPZERO_SV), pp_zero);
          fo[5] = CXCONDITIONAL_SV(sfomeg[0]*chi[ip], CXMAKE_SV_2ARG(fptype_sv(fptype(ip_zero)*sqm[sycl::abs(im_zero)]), FPZERO_SV), pp_zero);
      }
      else {
          const fptype_sv sqp0p3 = FPCONDITIONAL_SV(fptype(nsf)*FPSQRT(FPMAX(pvec0 + pvec3, 0.0)), FPZERO_SV, (pvec1 == 0.0) && (pvec2 == 0.0) && (pvec3 < 0.0));
          const cxtype_sv chi[2] = { CXMAKE_SV_2ARG( sqp0p3, FPZERO_SV ),
                                     CXCONDITIONAL_SV(CXMAKE_SV_2ARG(fptype(nh)*pvec1/sqp0p3, -pvec2/sqp0p3), CXMAKE_SV_2ARG(fptype(-nhel)*FPSQRT(2.0*pvec0), FPZERO_SV), sqp0p3 == 0.0)
                                   };
          if (nh == 1) {
              fo[2] = chi[0];
              fo[3] = chi[1];
              fo[4] = CXZERO_SV;
              fo[5] = CXZERO_SV;
          }
          else {
              fo[2] = CXZERO_SV;
              fo[3] = CXZERO_SV;
              fo[4] = chi[1];
              fo[5] = chi[0];
          }
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
  SYCL_EXTERNAL
  void opzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) {
      auto pvec3 = momenta[0].z;

      fo[0] = CXMAKE_SV_2ARG(fptype(nsf)*pvec3, fptype(nsf)*pvec3);
      fo[1] = CXZERO_SV;
      const cxtype_sv csqp0p3 = CXMAKE_SV_2ARG(fptype(nsf)*FPSQRT(2.0*pvec3), FPZERO_SV);
      fo[3] = CXZERO_SV;
      fo[4] = CXZERO_SV;
      if (nhel*nsf == 1) {
          fo[2] = csqp0p3;
          fo[5] = CXZERO_SV;
      }
      else {
          fo[2] = CXZERO_SV;
          fo[5] = csqp0p3;
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
  SYCL_EXTERNAL
  void omzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) {
      auto pvec3 = momenta[0].z;

      fo[0] = CXMAKE_SV_2ARG(fptype(-nsf)*pvec3, fptype(nsf)*pvec3); // remember pvec0 == -pvec3
      fo[1] = CXZERO_SV;
      const cxtype_sv chi1 = CXMAKE_SV_2ARG(fptype(-nhel)*FPSQRT(-2.0*pvec3), FPZERO_SV);
      if (nhel*nsf == 1) {
          fo[2] = CXZERO_SV;
          fo[3] = chi1;
          fo[4] = CXZERO_SV;
          fo[5] = CXZERO_SV;
      }
      else {
          fo[2] = CXZERO_SV;
          fo[3] = CXZERO_SV;
          fo[4] = chi1;
          fo[5] = CXZERO_SV;
      }
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction fo[6] from the input momenta[npar*4*nevt]
  // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
  SYCL_EXTERNAL
  void oxzxxx( const vector4* momenta,
               //const fptype fmass,           // ASSUME fermion mass==0
               const signed char nhel,         // input: -1 or +1 (helicity of fermion)
               const signed char nsf,          // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv fo[]
             ) {
      auto pvec0 = momenta[0].w;
      auto pvec1 = momenta[0].x;
      auto pvec2 = momenta[0].y;
      auto pvec3 = momenta[0].z;

      fo[0] = CXMAKE_SV_2ARG(fptype(nsf)*pvec0, fptype(nsf)*pvec3);
      fo[1] = CXMAKE_SV_2ARG(fptype(nsf)*pvec1, fptype(nsf)*pvec2);
      const fptype_sv sqp0p3 = fptype(nsf)*FPSQRT(pvec0 + pvec3);
      const cxtype_sv chi0 = CXMAKE_SV_2ARG(sqp0p3, FPZERO_SV);
      const cxtype_sv chi1 = CXMAKE_SV_2ARG(fptype(nhel*nsf)*pvec1/sqp0p3, -pvec2/sqp0p3);
      if (nhel*nsf == 1) {
          fo[2] = chi0;
          fo[3] = chi1;
          fo[4] = CXZERO_SV;
          fo[5] = CXZERO_SV;
      }
      else {
          fo[2] = CXZERO_SV;
          fo[3] = CXZERO_SV;
          fo[4] = chi1;
          fo[5] = chi0;
      }
  }

  //==========================================================================


  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
  SYCL_EXTERNAL INLINE
  void VVV1_0( const cxtype_sv V1[],
               const cxtype_sv V2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  SYCL_EXTERNAL INLINE
  void VVV1P0_1( const cxtype_sv V2[],
                 const cxtype_sv V3[],
                 const cxtype_sv COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  SYCL_EXTERNAL INLINE
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  SYCL_EXTERNAL INLINE
  void FFV1_1( const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               const fptype M1,
               const fptype W1,
               cxtype_sv F1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  SYCL_EXTERNAL INLINE
  void FFV1_2( const cxtype_sv F1[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               const fptype M2,
               const fptype W2,
               cxtype_sv F2[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  SYCL_EXTERNAL INLINE
  void FFV1P0_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype_sv COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV1_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV1P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV3_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV3P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV4_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex ) ALWAYS_INLINE;

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL INLINE
  void VVVV4P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] ) ALWAYS_INLINE;
! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Sep 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2023 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: O. Mattelaer (Mar 2020) for the MG5aMC CUDACPP plugin.
! Further modified by: O. Mattelaer, A. Valassi (2020-2023) for the MG5aMC CUDACPP plugin.
!==========================================================================

  //==========================================================================


  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6]
  SYCL_EXTERNAL
  void VVV1_0( const cxtype_sv V1[],
               const cxtype_sv V2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               cxtype_sv* vertex )
  {
    const fptype_sv P1[4] = { +CXREAL( V1[0] ), +CXREAL( V1[1] ), +CXIMAG( V1[1] ), +CXIMAG( V1[0] ) };
    const fptype_sv P2[4] = { +CXREAL( V2[0] ), +CXREAL( V2[1] ), +CXIMAG( V2[1] ), +CXIMAG( V2[0] ) };
    const fptype_sv P3[4] = { +CXREAL( V3[0] ), +CXREAL( V3[1] ), +CXIMAG( V3[1] ), +CXIMAG( V3[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP1 = ( V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5] );
    const cxtype_sv TMP2 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP3 = ( V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5] );
    const cxtype_sv TMP4 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP5 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv TMP7 = ( V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3] );
    const cxtype_sv TMP8 = ( V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3] );
    (*vertex) = COUP * ( TMP1 * ( - CXIMAGINARYI_SV * ( TMP0 ) + CXIMAGINARYI_SV * ( TMP2 ) ) + ( TMP3 * ( +CXIMAGINARYI_SV * ( TMP4 )- CXIMAGINARYI_SV * ( TMP5 ) ) + TMP6 * ( - CXIMAGINARYI_SV * ( TMP7 ) + CXIMAGINARYI_SV * ( TMP8 ) ) ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  SYCL_EXTERNAL
  void VVV1P0_1( const cxtype_sv V2[],
                 const cxtype_sv V3[],
                 const cxtype_sv COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype_sv V1[] )
  {
    const fptype_sv P2[4] = { +CXREAL( V2[0] ), +CXREAL( V2[1] ), +CXIMAG( V2[1] ), +CXIMAG( V2[0] ) };
    const fptype_sv P3[4] = { +CXREAL( V3[0] ), +CXREAL( V3[1] ), +CXIMAG( V3[1] ), +CXIMAG( V3[0] ) };
    V1[0] = + V2[0] + V3[0];
    V1[1] = + V2[1] + V3[1];
    const fptype_sv P1[4] = { -CXREAL( V1[0] ), -CXREAL( V1[1] ), -CXIMAG( V1[1] ), -CXIMAG( V1[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP2 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP4 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP5 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - CXIMAGINARYI_SV * W1 ) );
    V1[2] = denom * ( TMP6 * ( - CXIMAGINARYI_SV * ( P2[0] ) + CXIMAGINARYI_SV * ( P3[0] ) ) + ( V2[2] * ( - CXIMAGINARYI_SV * ( TMP0 ) + CXIMAGINARYI_SV * ( TMP2 ) ) + V3[2] * ( +CXIMAGINARYI_SV * ( TMP4 )- CXIMAGINARYI_SV * ( TMP5 ) ) ) );
    V1[3] = denom * ( TMP6 * ( - CXIMAGINARYI_SV * ( P2[1] ) + CXIMAGINARYI_SV * ( P3[1] ) ) + ( V2[3] * ( - CXIMAGINARYI_SV * ( TMP0 ) + CXIMAGINARYI_SV * ( TMP2 ) ) + V3[3] * ( +CXIMAGINARYI_SV * ( TMP4 )- CXIMAGINARYI_SV * ( TMP5 ) ) ) );
    V1[4] = denom * ( TMP6 * ( - CXIMAGINARYI_SV * ( P2[2] ) + CXIMAGINARYI_SV * ( P3[2] ) ) + ( V2[4] * ( - CXIMAGINARYI_SV * ( TMP0 ) + CXIMAGINARYI_SV * ( TMP2 ) ) + V3[4] * ( +CXIMAGINARYI_SV * ( TMP4 )- CXIMAGINARYI_SV * ( TMP5 ) ) ) );
    V1[5] = denom * ( TMP6 * ( - CXIMAGINARYI_SV * ( P2[3] ) + CXIMAGINARYI_SV * ( P3[3] ) ) + ( V2[5] * ( - CXIMAGINARYI_SV * ( TMP0 ) + CXIMAGINARYI_SV * ( TMP2 ) ) + V3[5] * ( +CXIMAGINARYI_SV * ( TMP4 )- CXIMAGINARYI_SV * ( TMP5 ) ) ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  SYCL_EXTERNAL
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               cxtype_sv* vertex )
  {
    const cxtype_sv TMP9 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) + ( F1[3] * ( F2[4] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    (*vertex) = COUP * - CXIMAGINARYI_SV * TMP9;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  SYCL_EXTERNAL
  void FFV1_1( const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               const fptype M1,
               const fptype W1,
               cxtype_sv F1[] )
  {
    F1[0] = + F2[0] + V3[0];
    F1[1] = + F2[1] + V3[1];
    const fptype_sv P1[4] = { -CXREAL( F1[0] ), -CXREAL( F1[1] ), -CXIMAG( F1[1] ), -CXIMAG( F1[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - CXIMAGINARYI_SV * W1 ) );
    F1[2] = denom * CXIMAGINARYI_SV * ( F2[2] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[2] * ( +CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[1] * (- one) * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +CXIMAGINARYI_SV * ( V3[2] + V3[5] ) ) + P1[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + M1 * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) );
    F1[3] = denom * (- CXIMAGINARYI_SV) * ( F2[2] * ( P1[0] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[1] * ( V3[2] - V3[5] ) + ( P1[2] * ( - CXIMAGINARYI_SV * ( V3[2] ) + CXIMAGINARYI_SV * ( V3[5] ) ) + P1[3] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * (- one) * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[2] * ( +CXIMAGINARYI_SV * ( V3[3] ) - V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + M1 * ( F2[4] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + F2[5] * ( -V3[2] + V3[5] ) ) ) );
    F1[4] = denom * (- CXIMAGINARYI_SV) * ( F2[4] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[2] * (- one) * ( +CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + ( F2[5] * ( P1[0] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[1] * ( -V3[2] + V3[5] ) + ( P1[2] * ( - CXIMAGINARYI_SV * ( V3[2] ) + CXIMAGINARYI_SV * ( V3[5] ) ) - P1[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + M1 * ( F2[2] * ( -V3[2] + V3[5] ) + F2[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) );
    F1[5] = denom * CXIMAGINARYI_SV * ( F2[4] * ( P1[0] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[1] * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +CXIMAGINARYI_SV * ( V3[2] + V3[5] ) ) + P1[3] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + ( F2[5] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P1[2] * ( - CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + M1 * ( F2[2] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  SYCL_EXTERNAL
  void FFV1_2( const cxtype_sv F1[],
               const cxtype_sv V3[],
               const cxtype_sv COUP,
               const fptype M2,
               const fptype W2,
               cxtype_sv F2[] )
  {
    F2[0] = + F1[0] + V3[0];
    F2[1] = + F1[1] + V3[1];
    const fptype_sv P2[4] = { -CXREAL( F2[0] ), -CXREAL( F2[1] ), -CXIMAG( F2[1] ), -CXIMAG( F2[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P2[0] * P2[0] ) - ( P2[1] * P2[1] ) - ( P2[2] * P2[2] ) - ( P2[3] * P2[3] ) - M2 * ( M2 - CXIMAGINARYI_SV * W2 ) );
    F2[2] = denom * CXIMAGINARYI_SV * ( F1[2] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * (- one) * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[3] ) - V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + ( F1[3] * ( P2[0] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[1] * ( -V3[2] + V3[5] ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[2] )- CXIMAGINARYI_SV * ( V3[5] ) ) + P2[3] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + M2 * ( F1[4] * ( V3[2] - V3[5] ) + F1[5] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) );
    F2[3] = denom * (- CXIMAGINARYI_SV) * ( F1[2] * ( P2[0] * (- one) * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[1] * ( V3[2] + V3[5] ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[2] + V3[5] ) ) - P2[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + ( F1[3] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + M2 * ( F1[4] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) - F1[5] * ( V3[2] + V3[5] ) ) ) );
    F2[4] = denom * (- CXIMAGINARYI_SV) * ( F1[4] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[2] * ( - CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[1] * (- one) * ( V3[2] + V3[5] ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[2] + V3[5] ) ) + P2[3] * ( V3[3]- CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + M2 * ( F1[2] * (- one) * ( V3[2] + V3[5] ) + F1[3] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) );
    F2[5] = denom * CXIMAGINARYI_SV * ( F1[4] * ( P2[0] * (- one) * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[1] * ( V3[2] - V3[5] ) + ( P2[2] * ( +CXIMAGINARYI_SV * ( V3[2] )- CXIMAGINARYI_SV * ( V3[5] ) ) + P2[3] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * ( -V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + ( P2[2] * (- one) * ( +CXIMAGINARYI_SV * ( V3[3] ) + V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + M2 * ( F1[2] * ( V3[3] + CXIMAGINARYI_SV * ( V3[4] ) ) + F1[3] * ( V3[2] - V3[5] ) ) ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  SYCL_EXTERNAL
  void FFV1P0_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype_sv COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )
  {
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -CXREAL( V3[0] ), -CXREAL( V3[1] ), -CXIMAG( V3[1] ), -CXIMAG( V3[0] ) };
    const cxtype_sv denom = COUP / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - CXIMAGINARYI_SV * W3 ) );
    V3[2] = denom * (- CXIMAGINARYI_SV) * ( F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3] );
    V3[3] = denom * (- CXIMAGINARYI_SV) * ( -F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] );
    V3[4] = denom * (- CXIMAGINARYI_SV) * ( - CXIMAGINARYI_SV * ( F1[2] * F2[5] + F1[5] * F2[2] ) + CXIMAGINARYI_SV * ( F1[3] * F2[4] + F1[4] * F2[3] ) );
    V3[5] = denom * (- CXIMAGINARYI_SV) * ( -F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2] );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV1_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex )
  {
    const cxtype_sv TMP10 = ( V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5] );
    const cxtype_sv TMP11 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP3 = ( V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    (*vertex) = COUP * ( - CXIMAGINARYI_SV * ( TMP6 * TMP10 ) + CXIMAGINARYI_SV * ( TMP3 * TMP11 ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV1P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] )
  {
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -CXREAL( V1[0] ), -CXREAL( V1[1] ), -CXIMAG( V1[1] ), -CXIMAG( V1[0] ) };
    const cxtype_sv TMP11 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - CXIMAGINARYI_SV * W1 ) );
    V1[2] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[2] ) + CXIMAGINARYI_SV * ( V3[2] * TMP11 ) );
    V1[3] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[3] ) + CXIMAGINARYI_SV * ( V3[3] * TMP11 ) );
    V1[4] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[4] ) + CXIMAGINARYI_SV * ( V3[4] * TMP11 ) );
    V1[5] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[5] ) + CXIMAGINARYI_SV * ( V3[5] * TMP11 ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV3_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex )
  {
    const cxtype_sv TMP1 = ( V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5] );
    const cxtype_sv TMP10 = ( V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5] );
    const cxtype_sv TMP12 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    (*vertex) = COUP * ( - CXIMAGINARYI_SV * ( TMP6 * TMP10 ) + CXIMAGINARYI_SV * ( TMP1 * TMP12 ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV3P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] )
  {
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -CXREAL( V1[0] ), -CXREAL( V1[1] ), -CXIMAG( V1[1] ), -CXIMAG( V1[0] ) };
    const cxtype_sv TMP12 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv TMP6 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - CXIMAGINARYI_SV * W1 ) );
    V1[2] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[2] ) + CXIMAGINARYI_SV * ( V2[2] * TMP12 ) );
    V1[3] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[3] ) + CXIMAGINARYI_SV * ( V2[3] * TMP12 ) );
    V1[4] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[4] ) + CXIMAGINARYI_SV * ( V2[4] * TMP12 ) );
    V1[5] = denom * ( - CXIMAGINARYI_SV * ( TMP6 * V4[5] ) + CXIMAGINARYI_SV * ( V2[5] * TMP12 ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV4_0( const cxtype_sv V1[],
                const cxtype_sv V2[],
                const cxtype_sv V3[],
                const cxtype_sv V4[],
                const cxtype_sv COUP,
                cxtype_sv* vertex )
  {
    const cxtype_sv TMP1 = ( V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5] );
    const cxtype_sv TMP11 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP12 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv TMP3 = ( V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5] );
    (*vertex) = COUP * ( - CXIMAGINARYI_SV * ( TMP3 * TMP11 ) + CXIMAGINARYI_SV * ( TMP1 * TMP12 ) );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6], V4[6]
  SYCL_EXTERNAL
  void VVVV4P0_1( const cxtype_sv V2[],
                  const cxtype_sv V3[],
                  const cxtype_sv V4[],
                  const cxtype_sv COUP,
                  const fptype M1,
                  const fptype W1,
                  cxtype_sv V1[] )
  {
    V1[0] = + V2[0] + V3[0] + V4[0];
    V1[1] = + V2[1] + V3[1] + V4[1];
    const fptype_sv P1[4] = { -CXREAL( V1[0] ), -CXREAL( V1[1] ), -CXIMAG( V1[1] ), -CXIMAG( V1[0] ) };
    const cxtype_sv TMP11 = ( V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5] );
    const cxtype_sv TMP12 = ( V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - CXIMAGINARYI_SV * W1 ) );
    V1[2] = denom * ( - CXIMAGINARYI_SV * ( V3[2] * TMP11 ) + CXIMAGINARYI_SV * ( V2[2] * TMP12 ) );
    V1[3] = denom * ( - CXIMAGINARYI_SV * ( V3[3] * TMP11 ) + CXIMAGINARYI_SV * ( V2[3] * TMP12 ) );
    V1[4] = denom * ( - CXIMAGINARYI_SV * ( V3[4] * TMP11 ) + CXIMAGINARYI_SV * ( V2[4] * TMP12 ) );
    V1[5] = denom * ( - CXIMAGINARYI_SV * ( V3[5] * TMP11 ) + CXIMAGINARYI_SV * ( V2[5] * TMP12 ) );
  }

  //--------------------------------------------------------------------------

} // end namespace MG5_sm

#endif // HelAmps_sm_H

