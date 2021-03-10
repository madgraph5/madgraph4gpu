# MG5aMC eemumu gridpack

This directory contains a (modified) MG5aMC gridpack for ee to mumu.

It contains all relevant Fortran code including MadEvent and all relevant data files.

## MG5aMC version and installation (2021.03.09)

The gridpack was created from a MG5aMC 2.9.2 installation in
  /eos/home-a/avalassi/2021/MG5aMC/2.9.2/

This is based on the latest 2.9.2 tarball:
-  https://launchpad.net/mg5amcnlo/2.0/2.9.x/+download/MG5_aMC_v2.9.2.tar.gz

It is complemented with the following patches for a few python3 bugs:
-  https://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/2.9.3/revision/307
-  https://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/2.9.3/revision/308
-  https://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/2.9.3/revision/310

## Python version

Python 3.8 was used from /cvmfs using the [setupPython38.sh](./setupPython38.sh) script.

## EEMUMU gridpack creation

The EEMUMU gridpack was created as follows

```
cd /eos/home-a/avalassi/2021/MG5aMC/2.9.2/
\rm -rf _EEMUMU
cat << EOF >> _EEMUMU.mg
generate e+ e- > mu+ mu-
output _EEMUMU
launch
set gridpack True
set ebeam1 750
set ebeam2 750
set ptl -1
set etal -1
set drll -1
{ time python3.8 ./bin/mg5_aMC -l DEBUG _EEMUMU.mg ; } >& _EEMUMU_log.txt
```

The [_EEMUMU_log.txt](_EEMUMU_log.txt) logfile has been copied to this directory.

The cross section computed during gridpack creation, originally displayed in _EEMUMU/HTML/run_01/results.html, is the following:

| Graph | Cross-Section (pb) | Error | Events (K) | Unwgt | Luminosity (pb^-1) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| _EEMUMU/SubProcesses/P1_ll_ll | 0.046903 | 0.000144 | 28.02 | 1592.0 | 0 |

The breakdown for the two MadEvent channels is the following:

| Graph | Cross-Section (pb) | Error | Events (K) | Unwgt | Luminosity (pb^-1) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [_EEMUMU/SubProcesses/P1_ll_ll/G1](_EEMUMU/run01_G1_log.txt) | 0.02333 | 0.000131 | 14.01 | 798.0 | 3.42e+04 |
| [_EEMUMU/SubProcesses/P1_ll_ll/G2](_EEMUMU/run01_G2_log.txt) | 0.02357 | 6.01e-05 | 14.01 | 794.0 | 3.37e+04 |

## EEMUMU gridpack unpacking

The EEMUMU gridpack was unpacked as follows in this directory, removing all libraries and binaries that can be rebuilt:
```
 ln -sf /eos/home-a/avalassi/2021/MG5aMC/2.9.2/ .
 mkdir eemumu
 cd eemumu
 cp -dpr ../2.9.2/_EEMUMU/run_01_gridpack.tar.gz .
 tar -xvzf run_01_gridpack.tar.gz 
 \rm run_01_gridpack.tar.gz 
 touch madevent/Events/.keepme
 \rm $(find . -name *.pyc)
 \rm madevent/Source/BIAS/dummy/*.o
 \rm madevent/lib/*.a
 \rm madevent/bin/madevent
 cd madevent/Source
 make clean
 cd -
```

## EEMUMU gridpack usage

First, rebuild 
```
 cd madevent/Source
 make
 cd ../SubProcesses/P1_ll_ll/
 make
 cd ../../..
```

Then generate 100 unweighted events
```
  ./run.sh 100 1
```

## EEMUMU gridpack: time performance for unweighted event generation (1)

On itscrd70, generating a variable number of unweighted events
```
  time ./run.sh <nunw> 1234
```
gives the following approximate performance (using gcc/gfortran 4.8.5, the default system compiler):

|   nunw |  real |  user |  sys |
|   ---: |  ---: |  ---: | ---: |
|   1000 |  2.1s |  0.9s | 0.2s |
|  10000 | 12.6s |  9.9s | 1.1s |
| 100000 | 90.9s | 79.4s | 7.6s |

Note that the 'real' time may fluctuate enormously from run to run with the same input parameters.

### Flamegraph (1)

A flamegraph produced with the following script for 100k events
```
  ../flgrAV time ./run.sh 100000 1234
```
is shown below. 

<img src="eemumu/eemumu_madevent_100Kunw.png"  width="1200"/>

The execution took 80s of user CPU time. The sampling rate is set at 1kHz, so 80k frames were collected. 
- Note that the `MATRIX1_` function, where matrix elements (ME) are computed, only took 4.8k frames, i.e. 4.8s of user CPU time.
- Note also that the unweighting efficiency, as seen in the gridpack creation step, is approximately 1.6k/28k. This means that, to generate 100k unweighted events, approximately 1.75M MEs were computed.
- This gives an approximate throughput of 1.75M/4.8s i.e. 3.7E5 MEs per second for this Fortran implementation. 

## EEMUMU gridpack: time performance for unweighted event generation (2)

Some performance counters and timers have later been added for the ME calculation, which is performed in the Fortran subroutine MATRIX1.

On itscrd70, generating a variable number of unweighted events
```
  time ./run.sh <nunw> 1234
```
gives the following approximate performance (using gcc/gfortran 4.8.5, the default system compiler):

|   nunw |   real |  user |  sys |    MATRIX1 calls |  MATRIX1 times | MATRIX1 throughputs |
|   ---: |   ---: |  ---: | ---: |             ---: |           ---: |                ---: |
|   1000 |   2.1s |  0.9s | 0.2s |       4019, 4019 | 0.015s, 0.015s |  2.67E5/s, 2.65E5/s |
|  10000 |  13.5s | 10.7s | 1.3s |   124019, 124019 | 0.471s, 0.459s |  2.63E5/s, 2.70E5/s |
| 100000 | 104.7s | 86.4s | 9.6s | 1020019, 1020019 | 3.776s, 3.792s |  2.70E5/s, 2.69E5/s |

Note that the 'real' time may fluctuate enormously from run to run with the same input parameters.

In the table above, each run yields two distinct sets of MATRIX1 metrics for the G1 and G2 subprocesses.
In fact, the madevent Fortran application is called twice, for the two different MadEvent channels.

Note that the timers have introduced a 10-20% overhead to user and sys CPU times, with respect to those in the previous table.

The number of matrix calls should indicate the number of phase space points before unweighting.
It is a bit puzzling that the ratio between unweighted and weighte events (i.e. the unweighting efficiency) seems to vary.
The MATRIX throughout however seems very stable, which seems to indicate that the numbers are correct.

Note, in particular, that the unweighting efficiency with 100k unweighted events is lower than previously estimated in the analysis of the flamegraph: 2.00M MEs (not 1.75M) were computed to generate 100k unweighted events. However, the ME throughput, 2.7E5, is lower than the 3.7E5 previously estimated, because the ME calculation took 7.6s instead of 4.8s, due to the overhead from the timers.

### Flamegraph (2)

A flamegraph produced with the following script for 100k events
```
  ../flgrAV time ./run.sh 100000 1234
```
is shown below. 

<img src="eemumu/eemumu_madevent_100Kunw2.png"  width="1200"/>

The execution took 86s of user CPU time. The sampling rate is set at 1kHz, so 86k frames were collected. 
- The `MATRIX1_` function, where matrix elements (ME) are computed, took 8.0k frames, i.e. 8.0s of user CPU time. Out of these, 3.0k frames, i.e. 3.0s, are an overhead from the counter start/stop functions.


