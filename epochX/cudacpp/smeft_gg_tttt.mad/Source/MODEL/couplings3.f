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
      GC_6(VECID) = -(MDL_COMPLEXI*G)
      GC_7(VECID) = G
      GC_8(VECID) = MDL_COMPLEXI*MDL_G__EXP__2
      END
