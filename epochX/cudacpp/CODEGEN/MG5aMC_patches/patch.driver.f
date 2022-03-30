diff --git a/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/driver.f b/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/driver.f
index 91e1f5b4..61a31a9f 100644
--- a/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/driver.f
+++ b/epochX/cudacpp/gg_tt.madonly/SubProcesses/P1_gg_ttx/driver.f
@@ -76,6 +76,7 @@ c      common/to_colstats/ncols,ncolflow,ncolalt,ic
 C-----
 C  BEGIN CODE
 C----- 
+      call counters_initialise()
       call cpu_time(t_before)
       CUMULATED_TIMING = t_before
 c
@@ -200,6 +201,7 @@ c      write(*,*) 'Final xsec: ',xsec
 
       rewind(lun)
 
+      call counters_finalise()
       close(lun)
       end
 
