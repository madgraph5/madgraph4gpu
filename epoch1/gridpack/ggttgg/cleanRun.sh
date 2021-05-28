#!/bin/bash 
# See "find . -newerct '5 minutes ago' -type f -printf '%CT %10s %p\n' | sort"

\rm -f ./events.lhe.gz
\rm -f ./madevent/Cards/grid_card.dat
\rm -rf ./madevent/Events/GridRun_*
\rm -f ./madevent/SubProcesses/P1_*/G*/results.dat
\rm -f ./madevent/SubProcesses/P1_*/G*/moffset.dat
\rm -f ./madevent/SubProcesses/P1_*/G*/input_sg.txt
\rm -f ./madevent/SubProcesses/P1_*/G*/counters_log.txt
\rm -f ./madevent/SubProcesses/P1_*/ajob1
\rm -f ./madevent/SubProcesses/P1_*/results.dat
\rm -f ./madevent/SubProcesses/randinit
\rm -f ./madevent/SubProcesses/results.dat
\rm -f ./madevent/index.html

# The following are only created if keeplog is true
\rm -f ./madevent/SubProcesses/P1_*/G*/ftn26
\rm -f ./madevent/SubProcesses/P1_*/G*/GridRun_*
