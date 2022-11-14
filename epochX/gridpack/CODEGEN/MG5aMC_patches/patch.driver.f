--- ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/driver.f
+++ ee_mumu.auto/madevent/SubProcesses/P1_ll_ll/driver.f.new
@@ -75,6 +75,7 @@
 C-----
 C  BEGIN CODE
 C----- 
+      call counters_initialise()
       call cpu_time(t_before)
       CUMULATED_TIMING = t_before
 c
@@ -194,6 +195,7 @@
       rewind(lun)
 
       close(lun)
+      call counters_finalise()
       end
 
 c     $B$ get_user_params $B$ ! tag for MadWeight
