diff --git a/epochX/cudacpp/gg_tt.mad/SubProcesses/makefile b/epochX/cudacpp/gg_tt.mad/SubProcesses/makefile
index cce95279..b8a2686d 100644
--- a/epochX/cudacpp/gg_tt.mad/SubProcesses/makefile
+++ b/epochX/cudacpp/gg_tt.mad/SubProcesses/makefile
@@ -1,6 +1,9 @@
 include ../../Source/make_opts
 FFLAGS+= -w
 
+# Enable the C preprocessor https://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
+FFLAGS+= -cpp 
+
 # Load additional dependencies of the bias module, if present
 ifeq (,$(wildcard ../bias_dependencies))
 BIASDEPENDENCIES =
@@ -46,8 +49,11 @@ SYMMETRY = symmetry.o idenparts.o
 
 # Binaries
 
-$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX)
-	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
+$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) counters.o
+	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp counters.o
+
+counters.o: counters.cpp timer.h
+	$(CXX) -std=c++11 -Wall -Wshadow -Wextra -c $< -o $@
 
 $(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
 	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
