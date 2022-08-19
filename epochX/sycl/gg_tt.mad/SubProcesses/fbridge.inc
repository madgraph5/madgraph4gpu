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
C Execute the matrix-element calculation "sequence" via a Bridge on SYCL device.
C - PBRIDGE: the memory address of the C++ Bridge
C - MOMENTA: the input 4-momenta Fortran array
C - GS:      the input Gs (running QCD coupling constant alphas)
C - MES:     the output matrix element Fortran array
C - CHANID:  the Feynman diagram to enhance in multi-channel mode if 1 to n (disable multi-channel if 0)
C
      INTERFACE
         SUBROUTINE FBRIDGESEQUENCE(PBRIDGE, MOMENTA, GS, MES, CHANID)
         INTEGER*8 PBRIDGE
         DOUBLE PRECISION MOMENTA(*)
         DOUBLE PRECISION GS(*)
         DOUBLE PRECISION MES(*)
         INTEGER*4 CHANID
         END SUBROUTINE FBRIDGESEQUENCE
      END INTERFACE
