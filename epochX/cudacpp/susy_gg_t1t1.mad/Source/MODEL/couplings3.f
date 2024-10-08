ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3( VECID)

      IMPLICIT NONE
      INTEGER VECID
      INCLUDE 'model_functions.inc'
      INCLUDE '../vector.inc'


      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_6(VECID) = -G
      GC_55(VECID) = -(MDL_COMPLEXI*G*MDL_I51X33)-MDL_COMPLEXI*G
     $ *MDL_I52X33
      GC_57(VECID) = -(MDL_COMPLEXI*G*MDL_I51X36)-MDL_COMPLEXI*G
     $ *MDL_I52X36
      GC_90(VECID) = MDL_COMPLEXI*MDL_G__EXP__2*MDL_I74X33
     $ +MDL_COMPLEXI*MDL_G__EXP__2*MDL_I75X33
      END
