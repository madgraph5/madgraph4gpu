
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

