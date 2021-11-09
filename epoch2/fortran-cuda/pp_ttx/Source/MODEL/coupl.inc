ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION G, ALL_G(17)
      COMMON/STRONG/ G, ALL_G

      DOUBLE COMPLEX GAL(2)
      COMMON/WEAK/ GAL

      DOUBLE PRECISION MU_R, ALL_MU_R(17)
      COMMON/RSCALE/ MU_R, ALL_MU_R

      DOUBLE PRECISION NF
      PARAMETER(NF=4D0)
      DOUBLE PRECISION NL
      PARAMETER(NL=2D0)

      DOUBLE PRECISION MDL_MT,MDL_MH,MDL_MZ,MDL_MW,MDL_MB,MDL_MTA

      COMMON/MASSES/ MDL_MT,MDL_MH,MDL_MZ,MDL_MW,MDL_MB,MDL_MTA


      DOUBLE PRECISION MDL_WW,MDL_WZ,MDL_WT,MDL_WH

      COMMON/WIDTHS/ MDL_WW,MDL_WZ,MDL_WT,MDL_WH


      DOUBLE COMPLEX GC_10(17), GC_11(17)

      COMMON/COUPLINGS/ GC_10, GC_11

