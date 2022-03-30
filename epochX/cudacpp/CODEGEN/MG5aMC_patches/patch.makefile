diff --git a/epochX/cudacpp/gg_tt.madonly/SubProcesses/makefile b/epochX/cudacpp/gg_tt.madonly/SubProcesses/makefile
index cce95279..b5bf594d 100644
--- a/epochX/cudacpp/gg_tt.madonly/SubProcesses/makefile
+++ b/epochX/cudacpp/gg_tt.madonly/SubProcesses/makefile
@@ -46,8 +46,11 @@ SYMMETRY = symmetry.o idenparts.o
 
 # Binaries
 
-$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX)
-	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
+$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) counters.o
+	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o
+
+counters.o: counters.cpp timer.h
+	$(CXX)  -std=c++11 -Wall -Wshadow -Wextra -c $<
 
 $(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
 	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
