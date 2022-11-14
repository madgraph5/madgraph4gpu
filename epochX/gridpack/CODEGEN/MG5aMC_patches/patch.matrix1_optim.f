--- ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f
+++ ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1_optim.f.new
@@ -90,7 +90,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
-
+      call counters_smatrix1_start()
       DO I=1,NEXTERNAL
         JC(I) = +1
       ENDDO
@@ -163,6 +163,7 @@
         ENDIF
       ENDIF
       ANS=ANS/DBLE(IDEN)
+      call counters_smatrix1_stop()
       END
 
 
@@ -253,6 +254,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       if (first) then
         first=.false.
         IF(ZERO.ne.0d0) fk_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -318,7 +320,7 @@
           enddo
         Enddo
       ENDDO
-
+      call counters_matrix1_stop()
       END
 
 
