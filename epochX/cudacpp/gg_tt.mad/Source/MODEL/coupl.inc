ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc, but this is not
C     included in coupl.inc because Fortran has no include guards
C     
C     NB: therefore coupl.inc introduces the definition of a separate
C     variable VECSIZE_MEMMAX_COUPL, whose value must be the same as
C     that of VECSIZE_MEMMAX (this is cross-checked in PROGRAM DRIVER)
C     

      INTEGER VECSIZE_MEMMAX_COUPL
      PARAMETER (VECSIZE_MEMMAX_COUPL=16384)

      DOUBLE PRECISION G, ALL_G(VECSIZE_MEMMAX_COUPL)
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R(VECSIZE_MEMMAX_COUPL)
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=4D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=2D0)

      DOUBLE PRECISION MDL_MB,MDL_MH,MDL_MT,MDL_MTA,MDL_MW,MDL_MZ

      COMMON/MASSES/ MDL_MB,MDL_MH,MDL_MT,MDL_MTA,MDL_MW,MDL_MZ


      DOUBLE PRECISION MDL_WH,MDL_WT,MDL_WW,MDL_WZ

      COMMON/WIDTHS/ MDL_WH,MDL_WT,MDL_WW,MDL_WZ


      DOUBLE COMPLEX GC_10(VECSIZE_MEMMAX_COUPL),
     $  GC_11(VECSIZE_MEMMAX_COUPL)

      COMMON/COUPLINGS/ GC_10, GC_11

