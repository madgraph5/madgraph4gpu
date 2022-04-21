diff --git a/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/matrix1.f b/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/matrix1.f
index 0a69788c..faf965bb 100644
--- a/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/matrix1.f
+++ b/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/matrix1.f
@@ -136,6 +136,7 @@ C     ----------
 C     BEGIN CODE
 C     ----------
 
+      call counters_smatrix1_start()
       NTRY(IMIRROR)=NTRY(IMIRROR)+1
       THIS_NTRY(IMIRROR) = THIS_NTRY(IMIRROR)+1
       DO I=1,NEXTERNAL
@@ -232,6 +233,7 @@ C       Include the Jacobian from helicity sampling
         WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
       ELSE
         ANS = 1D0
+        call counters_smatrix1_stop()
         RETURN
       ENDIF
       IF (ANS.NE.0D0.AND.(ISUM_HEL .NE. 1.OR.HEL_PICKED.EQ.-1)) THEN
@@ -276,6 +278,7 @@ C           Set right sign for ANS, based on sign of chosen helicity
         ENDIF
       ENDIF
       ANS=ANS/DBLE(IDEN)
+      call counters_smatrix1_stop()
       END
 
 
@@ -376,6 +379,7 @@ C     1 T(2,1,3,4)
 C     ----------
 C     BEGIN CODE
 C     ----------
+      call counters_matrix1_start()
       IF (FIRST) THEN
         FIRST=.FALSE.
         IF(ZERO.NE.0D0) FK_ZERO = SIGN(MAX(ABS(ZERO), ABS(ZERO
@@ -449,6 +453,7 @@ C     JAMPs contributing to orders ALL_ORDERS=1
         ENDDO
       ENDDO
 
+      call counters_matrix1_stop()
       END
 
       SUBROUTINE PRINT_ZERO_AMP_1()
