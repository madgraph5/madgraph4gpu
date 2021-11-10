diff --git a/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/driver.f b/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/driver.f
index 91e1f5b4..f34bd4b1 100644
--- a/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/driver.f
+++ b/epoch1/gridpack/eemumu/madevent/SubProcesses/P1_ll_ll/driver.f
@@ -73,6 +73,10 @@ c      common/to_colstats/ncols,ncolflow,ncolalt,ic
 
       include 'coupl.inc'
 
+c
+c     initialise C and C++ modules
+c
+      call counters_initialise()
 C-----
 C  BEGIN CODE
 C----- 
@@ -201,6 +205,10 @@ c      write(*,*) 'Final xsec: ',xsec
       rewind(lun)
 
       close(lun)
+c
+c     finalise C and C++ modules
+c
+      call counters_finalise()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
