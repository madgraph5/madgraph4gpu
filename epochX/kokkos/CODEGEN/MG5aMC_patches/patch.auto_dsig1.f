diff --git b/epochX/kokkos/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f a/epochX/kokkos/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
index 1734289bf..3c66b950f 100644
--- b/epochX/kokkos/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
+++ a/epochX/kokkos/gg_tt.mad/SubProcesses/P1_gg_ttx/auto_dsig1.f
@@ -76,13 +76,13 @@ C     Keep track of whether cuts already calculated for this event
 
       INTEGER SUBDIAG(MAXSPROC),IB(2)
       COMMON/TO_SUB_DIAG/SUBDIAG,IB
-      INCLUDE 'coupl.inc'
+      include 'vector.inc'
+      include 'coupl.inc'
       INCLUDE 'run.inc'
 C     Common blocks
       CHARACTER*7         PDLABEL,EPA_LABEL
       INTEGER       LHAID
       COMMON/TO_PDF/LHAID,PDLABEL,EPA_LABEL
-      INCLUDE 'vector.inc'
 C     jamp2 information
       DOUBLE PRECISION JAMP2(0:MAXFLOW, NB_PAGE_MAX)
       COMMON/TO_JAMPS/       JAMP2
@@ -219,7 +219,8 @@ C     ****************************************************
 C     
 C     CONSTANTS
 C     
-      INCLUDE '../../Source/vector.inc'
+      include 'vector.inc'
+      include 'coupl.inc'
       INCLUDE 'genps.inc'
       INCLUDE 'nexternal.inc'
       INCLUDE 'maxconfigs.inc'
@@ -288,7 +289,6 @@ C     jamp2 information
 
       INTEGER SUBDIAG(MAXSPROC),IB(2)
       COMMON/TO_SUB_DIAG/SUBDIAG,IB
-      INCLUDE 'coupl.inc'
       INCLUDE 'run.inc'
 
       DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE_MAX)
@@ -451,8 +451,10 @@ C
 
       USE OMP_LIB
 
+      IMPLICIT NONE
       INCLUDE 'nexternal.inc'
-      INCLUDE '../../Source/vector.inc'
+      include 'vector.inc'
+      include 'coupl.inc'
       INCLUDE 'maxamps.inc'
       DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE_MAX)
       DOUBLE PRECISION HEL_RAND(NB_PAGE_MAX)
@@ -462,22 +464,125 @@ C
       DOUBLE PRECISION JAMP2_MULTI(0:MAXFLOW, NB_PAGE_MAX)
 
       INTEGER IVEC
+      INTEGER IEXT
+
+      INTEGER                    ISUM_HEL
+      LOGICAL                    MULTI_CHANNEL
+      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL
+
+      LOGICAL FIRST_CHID
+      SAVE FIRST_CHID
+      DATA FIRST_CHID/.TRUE./
+      
+#ifdef MG5AMC_MEEXPORTER_KOKKOS
+      INCLUDE 'fbridge.inc'
+      INCLUDE 'fbridge_common.inc'
+      INCLUDE 'genps.inc'
+      INCLUDE 'run.inc'
+      DOUBLE PRECISION OUT2(NB_PAGE_MAX)
+      DOUBLE PRECISION CBYF1
+
+      INTEGER*4 NWARNINGS
+      SAVE NWARNINGS
+      DATA NWARNINGS/0/
+      
+      LOGICAL FIRST
+      SAVE FIRST
+      DATA FIRST/.TRUE./
+
+      IF( FBRIDGE_MODE .LE. 0 ) THEN ! (FortranOnly=0 or BothQuiet=-1 or BothDebug=-2)
+#endif
+        call counters_smatrix1multi_start( -1, nb_page_loop ) ! fortran=-1
+c!$OMP PARALLEL
+c!$OMP DO
+        DO IVEC=1, NB_PAGE_LOOP
+          CALL SMATRIX1(P_MULTI(0,1,IVEC),
+     &      hel_rand(IVEC),
+     &      channel,
+     &      out(IVEC),
+C    &      selected_hel(IVEC),
+     &      jamp2_multi(0,IVEC),
+     &      IVEC
+     &      )
+        ENDDO
+c!$OMP END DO
+c!$OMP END PARALLEL
+        call counters_smatrix1multi_stop( -1 ) ! fortran=-1
+#ifdef MG5AMC_MEEXPORTER_KOKKOS
+      ENDIF
 
+      IF( FBRIDGE_MODE .EQ. 1 .OR. FBRIDGE_MODE .LT. 0 ) THEN ! (CppOnly=1 or BothQuiet=-1 or BothDebug=-2)
+        IF( LIMHEL.NE.0 ) THEN
+          WRITE(6,*) 'ERROR! The kokkos bridge only supports LIMHEL=0'
+          STOP
+        ENDIF
+        IF ( FIRST ) THEN ! exclude first pass (helicity filtering) from timers (#461)
+          CALL FBRIDGESEQUENCE(FBRIDGE_PBRIDGE, P_MULTI, ALL_G, OUT2, 0) ! 0: multi channel disabled for helicity filtering
+          FIRST = .FALSE.
+c         ! This is a workaround for https://github.com/oliviermattelaer/mg5amc_test/issues/22 (see PR #486)
+          IF( FBRIDGE_MODE .EQ. 1 ) THEN ! (CppOnly=1 : SMATRIX1 is not called at all)
+            CALL RESET_CUMULATIVE_VARIABLE() ! mimic 'avoid bias of the initialization' within SMATRIX1
+          ENDIF
+        ENDIF
+        call counters_smatrix1multi_start( 0, nb_page_loop ) ! kokkos=0
+        IF ( .NOT. MULTI_CHANNEL ) THEN
+          CALL FBRIDGESEQUENCE(FBRIDGE_PBRIDGE,
+     &      P_MULTI, ALL_G, OUT2, 0) ! 0: multi channel disabled
+        ELSE
+          IF( SDE_STRAT.NE.1 ) THEN
+            WRITE(6,*) 'ERROR! The kokkos bridge requires SDE=1' ! multi channel single-diagram enhancement strategy
+            STOP
+          ENDIF
+          CALL FBRIDGESEQUENCE(FBRIDGE_PBRIDGE,
+     &      P_MULTI, ALL_G, OUT2, CHANNEL) ! 1-N: multi channel enabled
+        ENDIF
+        call counters_smatrix1multi_stop( 0 ) ! kokkos=0
+      ENDIF
+
+      IF( FBRIDGE_MODE .LT. 0 ) THEN ! (BothQuiet=-1 or BothDebug=-2)
+        DO IVEC=1, NB_PAGE_LOOP
+          CBYF1 = OUT2(IVEC)/OUT(IVEC) - 1
+          FBRIDGE_NCBYF1 = FBRIDGE_NCBYF1 + 1
+          FBRIDGE_CBYF1SUM = FBRIDGE_CBYF1SUM + CBYF1
+          FBRIDGE_CBYF1SUM2 = FBRIDGE_CBYF1SUM2 + CBYF1 * CBYF1
+          IF( CBYF1 .GT. FBRIDGE_CBYF1MAX ) FBRIDGE_CBYF1MAX = CBYF1
+          IF( CBYF1 .LT. FBRIDGE_CBYF1MIN ) FBRIDGE_CBYF1MIN = CBYF1
+          IF( FBRIDGE_MODE .EQ. -2 ) THEN ! (BothDebug=-2)
+            WRITE (*,'(I2,2E16.8,F23.11)')
+     &        IVEC, OUT(IVEC), OUT2(IVEC), 1+CBYF1
+          ENDIF
+          IF( ABS(CBYF1).GT.5E-5 .AND. NWARNINGS.LT.20 ) THEN
+            NWARNINGS = NWARNINGS + 1
+            WRITE (*,'(A,I2,A,I4,2E16.8,F23.11)')
+     &        'WARNING! (', NWARNINGS, '/20) Deviation more than 5E-5',
+     &        IVEC, OUT(IVEC), OUT2(IVEC), 1+CBYF1
+          ENDIF
+        END DO
+      ENDIF
+
+      IF( FBRIDGE_MODE .EQ. 1 .OR. FBRIDGE_MODE .LT. 0 ) THEN ! (CppOnly=1 or BothQuiet=-1 or BothDebug=-2)
+        DO IVEC=1, NB_PAGE_LOOP
+          OUT(IVEC) = OUT2(IVEC) ! use the kokkos ME instead of the fortran ME!
+        END DO
+      ENDIF
+
+      IF( FBRIDGE_MODE .EQ. 1 ) THEN ! (CppOnly=1 : SMATRIX1 is not called at all, JAMP2_MULTI is not filled)
+        DO IVEC=1, NB_PAGE_LOOP
+          JAMP2_MULTI(0,IVEC) = 2 ! workaround for https://github.com/oliviermattelaer/mg5amc_test/issues/14
+        END DO
+      ENDIF
+#endif
+
+      IF ( FIRST_CHID ) THEN
+        IF ( MULTI_CHANNEL ) THEN
+          WRITE (*,*) 'MULTI_CHANNEL = TRUE'
+        ELSE
+          WRITE (*,*) 'MULTI_CHANNEL = FALSE'
+        ENDIF
+        WRITE (*,*) 'CHANNEL_ID =', CHANNEL
+        FIRST_CHID = .FALSE.
+      ENDIF
 
-!$OMP PARALLEL
-!$OMP DO
-      DO IVEC=1, NB_PAGE_LOOP
-        CALL SMATRIX1(P_MULTI(0,1,IVEC),
-     &	                         hel_rand(IVEC),
-     &				 channel,
-     &				 out(IVEC),
-C       &				 selected_hel(IVEC),
-     &				 jamp2_multi(0,IVEC),
-     &				 IVEC
-     &				 )
-      ENDDO
-!$OMP END DO
-!$OMP END PARALLEL
       RETURN
       END
 
