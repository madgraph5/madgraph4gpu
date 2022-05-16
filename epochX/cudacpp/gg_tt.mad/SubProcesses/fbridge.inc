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
      
C
C Create a Bridge and return its pointer
C - PBRIDGE: the memory address of the C++ Bridge
C - NEVT:    the number of events in the Fortran arrays
C - NPAR:    the number of external particles in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY: remove it?)
C - NP4:     the number of momenta components, usually 4, in the Fortran arrays (KEPT FOR SANITY CHECKS ONLY: remove it?)
C
      INTERFACE
         SUBROUTINE FBRIDGECREATE(PBRIDGE, NEVT, NPAR, NP4)
         INTEGER*8 PBRIDGE
         INTEGER*4 NEVT
         INTEGER*4 NPAR
         INTEGER*4 NP4
         END SUBROUTINE FBRIDGECREATE
      END INTERFACE
      
C
C Delete a Bridge.
C - PBRIDGE: the memory address of the C++ Bridge
C
      INTERFACE
         SUBROUTINE FBRIDGEDELETE(PBRIDGE)
         INTEGER*8 PBRIDGE
         END SUBROUTINE FBRIDGEDELETE
      END INTERFACE
      
C
C Execute the matrix-element calculation "sequence" via a Bridge on GPU/CUDA or CUDA/C++.
C - PBRIDGE: the memory address of the C++ Bridge
C - MOMENTA: the input 4-momenta Fortran array
C - GS:      the input Gs (running QCD coupling constant alphas)
C - MES:     the output matrix element Fortran array
C
      INTERFACE
         SUBROUTINE FBRIDGESEQUENCE(PBRIDGE, MOMENTA, GS, MES)
         INTEGER*8 PBRIDGE
         DOUBLE PRECISION MOMENTA(*)
         DOUBLE PRECISION GS(*)
         DOUBLE PRECISION MES(*)
         END SUBROUTINE FBRIDGESEQUENCE
      END INTERFACE
