diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
index 018ed7a6..e40d30e2 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
@@ -281,7 +281,6 @@ C     jamp2 information
       DOUBLE PRECISION HEL_RAND(NB_PAGE)
       INTEGER SELECTED_HEL(NB_PAGE)
 
-
 C     
 C     local
 C     
@@ -443,10 +442,15 @@ C
       DOUBLE PRECISION JAMP2_MULTI(0:MAXFLOW, NB_PAGE)
 
       INTEGER IVEC
+      INTEGER IEXT
 
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      INCLUDE 'fbridge.inc'
+      DOUBLE PRECISION OUT2(NB_PAGE)
+#endif
 
-!$OMP PARALLEL
-!$OMP DO
+c!$OMP PARALLEL
+c!$OMP DO
       DO IVEC=1, NB_PAGE
         CALL SMATRIX1(P_MULTI(0,1,IVEC),
      &	                         hel_rand(IVEC),
@@ -457,8 +461,21 @@ C       &				 selected_hel(IVEC),
      &				 IVEC
      &				 )
       ENDDO
-!$OMP END DO
-!$OMP END PARALLEL
+c!$OMP END DO
+c!$OMP END PARALLEL
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGESEQUENCE(MEEXPORTER_PBRIDGE, P_MULTI, OUT2)
+      DO IVEC=1, NB_PAGE
+c       DO IEXT=1, NEXTERNAL
+c         WRITE (*,*) P_MULTI(0,IEXT,IVEC), P_MULTI(1,IEXT,IVEC),
+c    &                P_MULTI(2,IEXT,IVEC), P_MULTI(3,IEXT,IVEC)
+c       END DO
+        WRITE (*,*) IVEC, OUT(IVEC), OUT2(IVEC), OUT2(IVEC)/OUT(IVEC)
+c       WRITE (*,*)
+      END DO
+#endif
+
       RETURN
       END
 
