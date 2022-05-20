ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP()

      IMPLICIT NONE
      DOUBLE PRECISION PI, ZERO
      LOGICAL READLHA
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'model_functions.inc'
      LOGICAL UPDATELOOP
      COMMON /TO_UPDATELOOP/UPDATELOOP
      INCLUDE 'input.inc'
      include 'vector.inc'
      include 'coupl.inc' ! NB must also include vector.inc
      READLHA = .TRUE.
      INCLUDE 'intparam_definition.inc'
      CALL COUP1()
C     
couplings needed to be evaluated points by points
C     
      CALL COUP2(1)

      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM(VECID)

      IMPLICIT NONE
      INTEGER VECID
      DOUBLE PRECISION PI, ZERO
      LOGICAL READLHA
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      LOGICAL UPDATELOOP
      COMMON /TO_UPDATELOOP/UPDATELOOP
      INCLUDE 'model_functions.inc'
      INCLUDE 'input.inc'
      include 'vector.inc'
      include 'coupl.inc' ! NB must also include vector.inc
      READLHA = .FALSE.

      INCLUDE 'intparam_definition.inc'


C     
couplings needed to be evaluated points by points
C     
      ALL_G(VECID) = G
      CALL COUP2(VECID)

      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM2(MU_R2,AS2 ,VECID)

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      DOUBLE PRECISION MU_R2, AS2
      INTEGER VECID
      INCLUDE 'model_functions.inc'
      INCLUDE 'input.inc'
      include 'vector.inc'
      include 'coupl.inc' ! NB must also include vector.inc

      IF (MU_R2.GT.0D0) MU_R = MU_R2
      G = SQRT(4.0D0*PI*AS2)
      AS = AS2

      CALL UPDATE_AS_PARAM(VECID)


      RETURN
      END

