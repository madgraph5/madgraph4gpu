sample code generated with command sequence

```
add process p p > t t~
output pp_ttx --vector_size=16
launch
set fixed_scale 100
#set dynamical_scale_choice 3
set use_syst F
set mc_grouped_subproc F
set global_flag "-O3 -funroll-loops"
set aloha_flag "-O3 -ffast-math"
#set global_flag "-O1 -g -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall"
set nevents 10k
set etal 1
```
 
using MG5 branch 

```
bzr checkout lp:~maddevelopers/mg5amcnlo/3.1.1_lo_vectorization
```


