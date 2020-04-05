# Performance Data Files

## Files 
Each performance run contains 4 files, distinguished by _\<date>_ and _\<run#>_ 

* _\<date>-params-\<run#>.txt_: The parameters used as inputs for the performance run, each line contains three numbers, i.e. <#blocks> <#threads> <#iterations>
* _\<date>-perf-test-\<run#>.txt_: The output of the performance runs, usually done via "./check.ext -p \<#blocks> \<#threads> \<#iterations>
* _\<date>-nvidia-smi-dmon-\<run#>.txt_: The output of "nvidia-smi dmon ..." during the performance run
* _\<date>-nvidia-smi-pmon-\<run#>.txt_: The output of "nvidia-smi pmon ..." during the performance run

## Performance runs

sw version | date | run # | # events / execution | max threads | # configs | Comment
--- | --- | --- | --- | --- | --- | --- 
 | 20200402 | run1 | 6.3 * 10^6 | 384 | | 16 iterations fixed, try max number of threads / block
 | 20200402 | run2 | 6.3 * 10^6 | 384 | 108 | full mesh of configs 
 | 20200404 | run1 | 6.3 * 10^6 | 256 | 108 | step size of threads is power of 2 
 | 20200405 | run1 | 6.3 * 10^6 | 256 | 108 | re-run 20200404-*-run1, check variance
