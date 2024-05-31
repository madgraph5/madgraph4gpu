ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1( )

      IMPLICIT NONE

      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE '../vector.inc'
      INCLUDE 'coupl.inc'
      GC_38 = -(MDL_CPWWWL2*MDL_CW*MDL_EE__EXP__3*MDL_COMPLEXI)
     $ /(2.000000D+06*MDL_SW__EXP__3)
      GC_39 = (3.000000D+00*MDL_CW*MDL_CWWWL2*MDL_EE__EXP__3
     $ *MDL_COMPLEXI)/(2.000000D+06*MDL_SW__EXP__3)
      GC_109 = (MDL_CPWL2*MDL_CW*MDL_EE__EXP__3*MDL_COMPLEXI
     $ *MDL_VEV__EXP__2)/(4.000000D+06*MDL_SW__EXP__3)
      GC_110 = -(MDL_CW*MDL_CWL2*MDL_EE__EXP__3*MDL_COMPLEXI
     $ *MDL_VEV__EXP__2)/(8.000000D+06*MDL_SW__EXP__3)
      GC_113 = (MDL_CBL2*MDL_EE__EXP__3*MDL_COMPLEXI*MDL_VEV__EXP__2)
     $ /(8.000000D+06*MDL_CW*MDL_SW)
      GC_114 = (MDL_CPWL2*MDL_EE__EXP__3*MDL_COMPLEXI*MDL_VEV__EXP__2)
     $ /(4.000000D+06*MDL_CW*MDL_SW)
      GC_115 = -(MDL_CWL2*MDL_EE__EXP__3*MDL_COMPLEXI*MDL_VEV__EXP__2)
     $ /(8.000000D+06*MDL_CW*MDL_SW)
      GC_145 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM1X1)/(MDL_SW
     $ *MDL_SQRT__2)
      END
