# Performance Data Files

## Files 
Each performance run contains 4 files, distinguished by _\<date>_ and _\<run#>_ 

* _\<date>-params-\<run#>.txt_: The parameters used as inputs for the performance run, each line contains three numbers, i.e. <#blocks> <#threads> <#iterations>
* _\<date>-perf-test-\<run#>.txt_: The output of the performance runs, usually done via "./check.ext -p \<#blocks> \<#threads> \<#iterations>
* _\<date>-nvidia-smi-dmon-\<run#>.txt_: The output of "nvidia-smi dmon ..." during the performance run
* _\<date>-nvidia-smi-pmon-\<run#>.txt_: The output of "nvidia-smi pmon ..." during the performance run

## standalone_cpu 

Numbers below are are taken from a performance test for the standalone_cpu numbers of the e+e- -> mu+m- test for 6.3 Mio produced events (6,291,455 events to be exact). 

Test machine: pmpe12.cern.ch

 | value
--- | --- 
model | Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz
num cores | 32 (test execution is single threaded)
memory | 64 GB 

Test results (6 April 2020)

 | value
--- | --- 
TotalTimeInWaveFuncs | 1.630552e+01
MeanTimeinWaveFuncs | 6.479229e-07
StdDevWaveFuncs | 6.883308e-07
MinTimeInWaveFuncs | 6.050000e-07
MaxTimeInWaveFuncs | 1.832447e-03


## Performance runs (gpu)

sw version | date | run # | processor type | # events / execution | max threads | # configs | Comment
--- | --- | --- | --- | --- | --- | --- | ---
22-chi | 20200402 | run1 | GV100GL? | 6.3 * 10^6 | 384 | | 16 iterations fixed, try max number of threads / block
22-chi | 20200402 | run2 | GV100GL? | 6.3 * 10^6 | 384 | 108 | full mesh of configs 
22-chi | 20200404 | run1 | GV100GL? | 6.3 * 10^6 | 256 | 108 | step size of threads is power of 2 
22-chi | 20200405 | run1 | GV100GL | 6.3 * 10^6 | 256 | 108 | re-run 20200404-*-run1, check variance
