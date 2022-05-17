C
C Common blocks for internal communication between MadEvent components
C - MEEXPORTER_MODE: the operation mode of MEs from Fortran and from the C++ Bridge
C - MEEXPORTER_*CBYF*: statistics for the CudaCpp by Fortran ratios of MEs
C
      INTEGER*4 MEEXPORTER_MODE ! (CppOnly=1, FortranOnly=0, BothQuiet=-1, BothDebug=-2)
      COMMON/TO_MEEXPORTER_MODE/MEEXPORTER_MODE
      INTEGER*8 MEEXPORTER_NCBYF1 ! number of entries for ME ratio-1
      DOUBLE PRECISION MEEXPORTER_CBYF1SUM, MEEXPORTER_CBYF1SUM2, ! sum/sum2/min/max ME ratio-1
     &  MEEXPORTER_CBYF1MIN, MEEXPORTER_CBYF1MAX
      COMMON/TO_MEEXPORTER_CBYF1/MEEXPORTER_NCBYF1,
     &  MEEXPORTER_CBYF1SUM, MEEXPORTER_CBYF1SUM2,
     &  MEEXPORTER_CBYF1MIN, MEEXPORTER_CBYF1MAX

C
C Common block for internal communication between MadEvent components
C - MEEXPORTER_PBRIDGE: the memory address of the C++ Bridge
C
      INTEGER*8 MEEXPORTER_PBRIDGE ! 64bit memory address
      COMMON/TO_MEEXPORTER_PBRIDGE/MEEXPORTER_PBRIDGE
