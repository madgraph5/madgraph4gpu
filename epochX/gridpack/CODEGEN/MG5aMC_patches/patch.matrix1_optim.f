diff --git a/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f b/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f
index 1c149fda..f7613001 100644
--- a/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f
+++ b/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f
@@ -253,6 +253,7 @@ C     1 ColorOne()
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       if (first) then
         first=.false.
         IF(ZERO.ne.0d0) fk_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -317,7 +318,7 @@ C     ----------
           enddo
         Enddo
       ENDDO
-
+      call counters_matrix1_stop()
       END
 
 
