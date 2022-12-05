
These process directories were generated using this version of MadGraph:

bzr branch lp:~maddevelopers/mg5amcnlo/2.7.0_gpu

Running this command:
```bash
generate g g > t t~ g g g
output standalone gg_ttggg
```

module load cmake/3.22.1 gcc/9.2.0

you can run the `run_test.py` file to execute all processes in order. 

I made changes to the `check_sa.f` to run the matrix element calculation 1000 times, calculate and print the average and stddev of the run time.

When I run this on a SkyLake (8180) with `mpirun -n 56 ...` (so that I fill up all physical cores) I see these numbers from the `run_test.py` file:
```
ee_mumu/SubProcesses/P1_epem_mupmum: 2.624715e-06 +/- 1.736501e-07  rate:  2.133565e+07
       gg_tt/SubProcesses/P1_gg_ttx: 1.609207e-05 +/- 3.800823e-07  rate:  3.479974e+06
     gg_ttg/SubProcesses/P1_gg_ttxg: 1.230965e-04 +/- 2.861196e-06  rate:  4.549275e+05
   gg_ttgg/SubProcesses/P1_gg_ttxgg: 1.444539e-03 +/- 2.689728e-05  rate:  3.876670e+04
 gg_ttggg/SubProcesses/P1_gg_ttxggg: 3.449992e-02 +/- 2.272253e-04  rate:  1.623192e+03
```