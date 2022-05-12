diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
index 91e1f5b4..4cfdb381 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
@@ -73,11 +73,23 @@ c      common/to_colstats/ncols,ncolflow,ncolalt,ic
 
       include 'coupl.inc'
 
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      INCLUDE '../../Source/vector.inc'
+      INCLUDE 'fbridge.inc'
+#endif
 C-----
 C  BEGIN CODE
 C----- 
       call cpu_time(t_before)
       CUMULATED_TIMING = t_before
+
+      CALL COUNTERS_INITIALISE()
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGECREATE(MEEXPORTER_PBRIDGE, NB_PAGE, NEXTERNAL, 4) ! this must be at the beginning as it initialises the CUDA device
+      MEEXPORTER_MODE = -1 ! (CppOnly=1, FortranOnly=0, BothQuiet=-1, BothDebug=-2)
+      MEEXPORTER_CBYFMAX = 0
+      MEEXPORTER_CBYFMIN = 1D100
+#endif
 c
 c     Read process number
 c
@@ -199,8 +211,18 @@ c      call sample_result(xsec,xerr)
 c      write(*,*) 'Final xsec: ',xsec
 
       rewind(lun)
-
       close(lun)
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGEDELETE(MEEXPORTER_PBRIDGE) ! this must be at the end as it shuts down the CUDA device
+      IF( MEEXPORTER_MODE .LE. -1 ) THEN ! (BothQuiet=-1 or BothDebug=-2)
+        WRITE(*,*) 'ME ratio CudaCpp/Fortran: MIN = ', MEEXPORTER_CBYFMIN
+        WRITE(*,*) 'ME ratio CudaCpp/Fortran: MAX = ', MEEXPORTER_CBYFMAX
+        WRITE(*,*) 'ME ratio CudaCpp/Fortran: 1-MIN = ', 1-MEEXPORTER_CBYFMIN
+        WRITE(*,*) 'ME ratio CudaCpp/Fortran: MAX-1 = ', MEEXPORTER_CBYFMAX-1
+      ENDIF
+#endif
+      CALL COUNTERS_FINALISE()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
