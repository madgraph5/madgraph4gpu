//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

namespace MG5_sm
{
  //--------------------------------------------------------------------------

  __device__
  void ixxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void ipzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void imzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void ixzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void vxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,
               const int nsv,
               cxtype* vc,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void sxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype,             // WARNING: "smass" unused???
               const int,                // WARNING: "nhel" unused???
               const int nss,
               cxtype sc[3],             // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar
  //--------------------------------------------------------------------------

  __device__
  void oxxxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void opzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void omzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar
  //--------------------------------------------------------------------------

  __device__
  void oxzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar );         // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__
  void FFV1_0( const cxtype F1[],    // input: wavefunction1[6]
               const cxtype F2[],    // input: wavefunction2[6]
               const cxtype V3[],    // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex );     // output: amplitude
  
  //--------------------------------------------------------------------------

  __device__
  void FFV1P0_3( const cxtype F1[],   // input: wavefunction1[6]
                 const cxtype F2[],   // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] );       // output: wavefunction3[6]
  
  //--------------------------------------------------------------------------

  __device__
  void FFV2_0( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype V3[],   // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex );    // output: amplitude

  //--------------------------------------------------------------------------

  __device__
  void FFV2_3( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] );       // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__
  void FFV4_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype* vertex );
  
  //--------------------------------------------------------------------------

  __device__
  void FFV4_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] );

  //--------------------------------------------------------------------------

  __device__
  void FFV2_4_0( const cxtype F1[],    // input: wavefunction1[6]
                 const cxtype F2[],    // input: wavefunction2[6]
                 const cxtype V3[],    // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype* vertex );     // output: amplitude
  
  //--------------------------------------------------------------------------

  __device__
  void FFV2_4_3( const cxtype F1[],   // input: wavefunction1[6]
                 const cxtype F2[],   // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] );       // output: wavefunction3[6]

  //--------------------------------------------------------------------------

} // end namespace

#endif // HelAmps_sm_H
