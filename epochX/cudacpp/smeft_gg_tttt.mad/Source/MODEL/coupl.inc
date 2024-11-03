ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     
C     NB: VECSIZE_MEMMAX is defined in vector.inc
C     NB: vector.inc must be included before coupl.inc
C     

      DOUBLE PRECISION G, ALL_G(VECSIZE_MEMMAX)
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R(VECSIZE_MEMMAX)
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=0D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=0D0)

      DOUBLE PRECISION MDL_MB,MDL_MC,MDL_MD,MDL_MH,MDL_MH1,MDL_MMU
     $ ,MDL_MS,MDL_MT,MDL_MT1,MDL_MTA,MDL_MU,MDL_MW,MDL_MW1,MDL_MZ
     $ ,MDL_MZ1,MDL_ME

      COMMON/MASSES/ MDL_MB,MDL_MC,MDL_MD,MDL_MH,MDL_MH1,MDL_MMU
     $ ,MDL_MS,MDL_MT,MDL_MT1,MDL_MTA,MDL_MU,MDL_MW,MDL_MW1,MDL_MZ
     $ ,MDL_MZ1,MDL_ME


      DOUBLE PRECISION MDL_WH,MDL_WH1,MDL_WT,MDL_WT1,MDL_WW,MDL_WW1
     $ ,MDL_WZ,MDL_WZ1

      COMMON/WIDTHS/ MDL_WH,MDL_WH1,MDL_WT,MDL_WT1,MDL_WW,MDL_WW1
     $ ,MDL_WZ,MDL_WZ1


      DOUBLE COMPLEX GC_6(VECSIZE_MEMMAX), GC_7(VECSIZE_MEMMAX),
     $  GC_8(VECSIZE_MEMMAX)

      COMMON/COUPLINGS/ GC_6, GC_7, GC_8

