diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
index 91e1f5b4..9daa40c8 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
@@ -73,11 +73,20 @@ c      common/to_colstats/ncols,ncolflow,ncolalt,ic
 
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
+#endif
 c
 c     Read process number
 c
@@ -199,8 +208,12 @@ c      call sample_result(xsec,xerr)
 c      write(*,*) 'Final xsec: ',xsec
 
       rewind(lun)
-
       close(lun)
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGEDELETE(MEEXPORTER_PBRIDGE) ! this must be at the end as it shuts down the CUDA device
+#endif
+      CALL COUNTERS_FINALISE()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
