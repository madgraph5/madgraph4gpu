--- ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f
+++ ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f.new
@@ -253,6 +253,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       if (first) then
         first=.false.
         IF(ZERO.ne.0d0) fk_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -318,7 +319,7 @@
           enddo
         Enddo
       ENDDO
-
+      call counters_matrix1_stop()
       END
 
 

