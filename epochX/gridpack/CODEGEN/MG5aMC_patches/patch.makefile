diff --git a/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/makefile b/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/makefile
index 1fe4c746..8d16ce86 100644
--- a/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/makefile
+++ b/epochX/gridpack/ee_mumu.auto/madevent/SubProcesses/makefile
@@ -46,8 +46,11 @@ SYMMETRY = symmetry.o idenparts.o
 
 # Binaries
 
-$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX)
-	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES)
+$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) counters.o
+	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) counters.o
+
+counters.o: counters.cpp timer.h
+	$(CXX)  -std=c++11 -Wall -Wshadow -Wextra -c $<
 
 $(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
 	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES)
