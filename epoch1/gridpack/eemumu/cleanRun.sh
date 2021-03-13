#!/bin/bash 
# See "find . -newerct '5 minutes ago' -type f -printf '%CT %10s %p\n' | sort"

\rm -f ./events.lhe.gz
\rm -f ./madevent/Cards/grid_card.dat
\rm -rf ./madevent/Events/GridRun_*
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/results.dat
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/moffset.dat
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/input_sg.txt
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/counters_log.txt
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/results.dat
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/moffset.dat
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/input_sg.txt
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/counters_log.txt
\rm -f ./madevent/SubProcesses/P1_ll_ll/ajob1
\rm -f ./madevent/SubProcesses/P1_ll_ll/results.dat
\rm -f ./madevent/SubProcesses/randinit
\rm -f ./madevent/SubProcesses/results.dat
\rm -f ./madevent/index.html

# The following are only created if keeplog is true
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/ftn26
\rm -f ./madevent/SubProcesses/P1_ll_ll/G1/GridRun_*
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/ftn26
\rm -f ./madevent/SubProcesses/P1_ll_ll/G2/GridRun_*
