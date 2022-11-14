--- makefile.orig	2022-06-09 18:23:36.677775186 +0200
+++ makefile	2022-06-09 19:02:31.736865066 +0200
@@ -52,8 +52,11 @@
 
 # Binaries
 
-$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX)
-	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
+$(PROG): $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX) counters.o
+	$(FC) -o $(PROG) $(PROCESS) $(MATRIX) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) counters.o -fopenmp
+
+counters.o: counters.cpp timer.h
+	$(CXX)  -std=c++11 -Wall -Wshadow -Wextra -c $<
 
 $(PROG)_forhel: $(PROCESS) auto_dsig.o $(LIBS) $(MATRIX_HEL)
 	$(FC) -o $(PROG)_forhel $(PROCESS) $(MATRIX_HEL) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) -fopenmp
