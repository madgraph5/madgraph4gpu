diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
index 3aa1040b..f69ef9a5 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
@@ -454,10 +454,16 @@ C
       DOUBLE PRECISION JAMP2_MULTI(0:MAXFLOW, NB_PAGE)
 
       INTEGER IVEC
+      INTEGER IEXT
 
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      INCLUDE 'coupl.inc'
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
@@ -468,8 +474,21 @@ C       &				 selected_hel(IVEC),
      &				 IVEC
      &				 )
       ENDDO
-!$OMP END DO
-!$OMP END PARALLEL
+c!$OMP END DO
+c!$OMP END PARALLEL
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGESEQUENCE(MEEXPORTER_PBRIDGE, P_MULTI, ALL_G, OUT2)
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
 
