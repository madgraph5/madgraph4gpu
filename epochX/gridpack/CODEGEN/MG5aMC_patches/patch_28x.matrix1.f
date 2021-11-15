--- /data/avalassi/GPU2020/madgraph4gpuX/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1.f.new	2021-11-15 21:13:03.459641279 +0100
+++ /data/avalassi/GPU2020/madgraph4gpuX/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/matrix1.f	2021-11-15 21:13:42.689590640 +0100
@@ -123,7 +123,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
-
+      call counters_smatrix1_start()
       NTRY(IMIRROR)=NTRY(IMIRROR)+1
       THIS_NTRY(IMIRROR) = THIS_NTRY(IMIRROR)+1
       DO I=1,NEXTERNAL
@@ -261,6 +261,7 @@
         ENDIF
       ENDIF
       ANS=ANS/DBLE(IDEN)
+      call counters_smatrix1_stop()
       END
 
 
@@ -342,6 +343,7 @@
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       IF (FIRST) THEN
         FIRST=.FALSE.
         IF(ZERO.NE.0D0) FK_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -390,7 +392,7 @@
           ENDDO
         ENDDO
       ENDDO
-
+      call counters_matrix1_stop()
       END
 
 C     Set of functions to handle the array indices of the split orders
