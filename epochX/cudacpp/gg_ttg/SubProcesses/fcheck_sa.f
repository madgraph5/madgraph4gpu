      PROGRAM FCHECK_SA
      IMPLICIT NONE
      INCLUDE 'fsampler.inc'
      INCLUDE 'fbridge.inc'
      INTEGER*8 SAMPLER, BRIDGE ! 64bit memory addresses
      INTEGER NEVTMAX, NEXTERNAL, NP4
      PARAMETER(NEVTMAX=2048*256, NEXTERNAL=4, NP4=4)
      CHARACTER*4 ARG0, ARG1, ARG2, ARG3
      INTEGER NARG1, NARG2, NARG3
      INTEGER NEVT, NITER
      INTEGER IEVT, IITER
      DOUBLE PRECISION MOMENTA(0:NP4-1, NEXTERNAL, NEVTMAX) ! c-array momenta[nevt][nexternal][np4]
      DOUBLE PRECISION MES(NEVTMAX)
      DOUBLE PRECISION MES_SUM ! use REAL*16 for quadruple precision
C
C READ COMMAND LINE ARGUMENTS
C (NB: most errors will crash the program !)
C
      IF ( COMMAND_ARGUMENT_COUNT() == 3 ) THEN
        CALL GET_COMMAND_ARGUMENT(1,ARG1)
        CALL GET_COMMAND_ARGUMENT(2,ARG2)
        CALL GET_COMMAND_ARGUMENT(3,ARG3)
        READ (ARG1,'(I4)') NARG1
        READ (ARG2,'(I4)') NARG2
        READ (ARG3,'(I4)') NARG3
        WRITE(6,*) "GPUBLOCKS=  ", NARG1
        WRITE(6,*) "GPUTHREADS= ", NARG2
        WRITE(6,*) "NITERATIONS=", NARG3
        NEVT = NARG1 * NARG2
        NITER = NARG3
        IF ( NEVT > NEVTMAX ) THEN
          WRITE(6,*) "ERROR! NEVT>NEVTMAX"
          STOP
        ENDIF
      ELSE
        CALL GET_COMMAND_ARGUMENT(0,ARG0)
        WRITE(6,*) "Usage: ", TRIM(ARG0),
     &    " gpublocks gputhreads niterations"
        STOP
      ENDIF
C
C USE SAMPLER AND BRIDGE
C
      MES_SUM = 0
      CALL FBRIDGECREATE(BRIDGE, NEVT, NEXTERNAL, NP4) ! this must be at the beginning as it initialises the CUDA device
      CALL FSAMPLERCREATE(SAMPLER, NEVT, NEXTERNAL, NP4)
      DO IITER = 1, NITER
        CALL FSAMPLERSEQUENCE(SAMPLER, MOMENTA)
        CALL FBRIDGESEQUENCE(BRIDGE, MOMENTA, MES)
        DO IEVT = 1, NEVT
c         WRITE(6,*) 'MES', IEVT, MES(IEVT)
          MES_SUM = MES_SUM + MES(IEVT)
        END DO
c       WRITE(6,*)
      END DO
      CALL FSAMPLERDELETE(SAMPLER)
      CALL FBRIDGEDELETE(BRIDGE) ! this must be at the end as it shuts down the CUDA device
      WRITE(6,*) 'Average Matrix Element:', MES_SUM/NEVT/NITER
      END PROGRAM FCHECK_SA
