diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
index 018ed7a6..0accde15 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
@@ -281,7 +281,6 @@ C     jamp2 information
       DOUBLE PRECISION HEL_RAND(NB_PAGE)
       INTEGER SELECTED_HEL(NB_PAGE)
 
-
 C     
 C     local
 C     
@@ -443,7 +442,12 @@ C
       DOUBLE PRECISION JAMP2_MULTI(0:MAXFLOW, NB_PAGE)
 
       INTEGER IVEC
+      INTEGER IEXT
 
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      INCLUDE 'fbridge.inc'
+      DOUBLE PRECISION OUT2(NB_PAGE)
+#endif
 
 !$OMP PARALLEL
 !$OMP DO
@@ -459,6 +463,19 @@ C       &				 selected_hel(IVEC),
       ENDDO
 !$OMP END DO
 !$OMP END PARALLEL
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGESEQUENCE(MEEXPORTER_PBRIDGE, P_MULTI, OUT2)
+      DO IVEC=1, NB_PAGE
+        DO IEXT=1, NEXTERNAL
+          WRITE (*,*) P_MULTI(0,IEXT,IVEC), P_MULTI(1,IEXT,IVEC),
+     &                P_MULTI(2,IEXT,IVEC), P_MULTI(3,IEXT,IVEC)
+        END DO
+        WRITE (*,*) IVEC, OUT(IVEC), OUT2(IVEC), OUT2(IVEC)/OUT(IVEC)
+        WRITE (*,*)
+      END DO
+#endif
+
       RETURN
       END
 
