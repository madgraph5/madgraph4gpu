--- ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1.f
+++ ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1.f.new
@@ -342,6 +342,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       IF (FIRST) THEN
         FIRST=.FALSE.
         IF(ZERO.NE.0D0) FK_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -390,7 +391,7 @@
           ENDDO
         ENDDO
       ENDDO
-
+      call counters_matrix1_stop()
       END
 
 C     Set of functions to handle the array indices of the split orders
