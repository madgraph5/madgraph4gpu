--- ee_mumu.auto/madevent/SubProcesses/makefile
+++ ee_mumu.auto/madevent/SubProcesses/makefile.new
@@ -40,8 +40,11 @@
 
 # Binaries
 
-$(PROG): $(PROCESS) auto_dsig.o $(LIBS)
-	$(FC) -o $(PROG) $(PROCESS) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES)
+$(PROG): $(PROCESS) auto_dsig.o $(LIBS) counters.o
+	$(FC) -o $(PROG) $(PROCESS) $(LINKLIBS) $(LDFLAGS) $(BIASDEPENDENCIES) counters.o
+
+counters.o: counters.cpp timer.h
+	$(CXX)  -std=c++11 -Wall -Wshadow -Wextra -c $<
 
 gensym: $(SYMMETRY) configs.inc $(LIBDIR)libmodel.$(libext) $(LIBDIR)libgeneric.$(libext)
 	$(FC) -o gensym $(SYMMETRY) -L../../lib/ -lmodel -lgeneric -lpdf $(llhapdf) $(LDFLAGS)
