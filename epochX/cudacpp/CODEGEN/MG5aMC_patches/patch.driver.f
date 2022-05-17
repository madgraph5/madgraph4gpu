diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
index 91e1f5b4..5648f6b0 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/P1_gg_ttx/driver.f
@@ -73,11 +73,27 @@ c      common/to_colstats/ncols,ncolflow,ncolalt,ic
 
       include 'coupl.inc'
 
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      INCLUDE '../../Source/vector.inc'
+      INCLUDE 'fbridge.inc'
+      INCLUDE 'fbridge_common.inc'
+#endif
 C-----
 C  BEGIN CODE
 C----- 
       call cpu_time(t_before)
       CUMULATED_TIMING = t_before
+
+      CALL COUNTERS_INITIALISE()
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGECREATE(FBRIDGE_PBRIDGE, NB_PAGE, NEXTERNAL, 4) ! this must be at the beginning as it initialises the CUDA device
+      FBRIDGE_MODE = -1 ! (CppOnly=1, FortranOnly=0, BothQuiet=-1, BothDebug=-2)
+      FBRIDGE_NCBYF1 = 0
+      FBRIDGE_CBYF1SUM = 0
+      FBRIDGE_CBYF1SUM2 = 0
+      FBRIDGE_CBYF1MAX = -1D100
+      FBRIDGE_CBYF1MIN = 1D100
+#endif
 c
 c     Read process number
 c
@@ -132,7 +148,8 @@ c   If CKKW-type matching, read IS Sudakov grid
           exit
  30       issgridfile='../'//issgridfile
           if(i.eq.5)then
-            print *,'ERROR: No Sudakov grid file found in lib with ickkw=2'
+            print *,
+     &        'ERROR: No Sudakov grid file found in lib with ickkw=2'
             stop
           endif
         enddo
@@ -199,8 +216,33 @@ c      call sample_result(xsec,xerr)
 c      write(*,*) 'Final xsec: ',xsec
 
       rewind(lun)
-
       close(lun)
+
+#ifdef MG5AMC_MEEXPORTER_CUDACPP
+      CALL FBRIDGEDELETE(FBRIDGE_PBRIDGE) ! this must be at the end as it shuts down the CUDA device
+      IF( FBRIDGE_MODE .LE. -1 ) THEN ! (BothQuiet=-1 or BothDebug=-2)
+        WRITE(*,'(a,f10.8,a,e8.2)')
+     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: MIN = ',
+     &    FBRIDGE_CBYF1MIN + 1, ' = 1 - ', -FBRIDGE_CBYF1MIN
+        WRITE(*,'(a,f10.8,a,e8.2)')
+     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: MAX = ',
+     &    FBRIDGE_CBYF1MAX + 1, ' = 1 + ', FBRIDGE_CBYF1MAX
+        WRITE(*,'(a,i6)')
+     &    ' [MERATIOS] ME ratio CudaCpp/Fortran: NENTRIES = ',
+     &    FBRIDGE_NCBYF1
+c        WRITE(*,'(a,e8.2)')
+c    &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: AVG = ',
+c    &    FBRIDGE_CBYF1SUM / FBRIDGE_NCBYF1
+c       WRITE(*,'(a,e8.2)')
+c    &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: STD = ',
+c    &    SQRT( FBRIDGE_CBYF1SUM2 / FBRIDGE_NCBYF1 ) ! ~standard deviation
+        WRITE(*,'(a,e8.2,a,e8.2)')
+     &    ' [MERATIOS] ME ratio CudaCpp/Fortran - 1: AVG = ',
+     &    FBRIDGE_CBYF1SUM / FBRIDGE_NCBYF1, ' +- ',
+     &    SQRT( FBRIDGE_CBYF1SUM2 ) / FBRIDGE_NCBYF1 ! ~standard error
+      ENDIF
+#endif
+      CALL COUNTERS_FINALISE()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
