diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/reweight.f b/epochX/cudacpp/gg_tt.mad/SubProcesses/reweight.f
index e8dc44c8..daaab251 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/reweight.f
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/reweight.f
@@ -1793,7 +1793,7 @@ C
       include 'run.inc'
       include 'nexternal.inc'
       include 'coupl.inc'
-#      include 'maxparticles.inc'
+C      include 'maxparticles.inc'
       
       double precision all_p(4*maxdim/3+14,*), all_wgt(*)
       double precision all_q2fact(2,*)
