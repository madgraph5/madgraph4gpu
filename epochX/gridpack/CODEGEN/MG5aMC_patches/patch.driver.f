--- ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/driver.f
+++ ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/driver.f.new
@@ -73,6 +73,10 @@
 
       include 'coupl.inc'
 
+c
+c     initialise C and C++ modules
+c
+      call counters_initialise()
 C-----
 C  BEGIN CODE
 C----- 
@@ -201,6 +205,10 @@
       rewind(lun)
 
       close(lun)
+c
+c     finalise C and C++ modules
+c
+      call counters_finalise()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
