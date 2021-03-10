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

|   nunw |   real |  user |  sys |    MATRIX1 calls |
|   ---: |   ---: |  ---: | ---: |             ---: |
|   1000 |   1.9s |  0.9s | 0.2s |       4019, 4019 |
|  10000 |  12.6s |  9.8s | 1.1s |   124019, 124019 |
| 100000 |  90.9s | 79.4s | 7.6s | 1020019, 1020019 |

Note that the 'real' time may fluctuate enormously from run to run with the same input parameters.

The MATRIX1 performance counters indicate how many times the script triggers the Fortran subroutine MATRIX1, where matrix elements (MEs) are computed. These counters, which have been added in epoch1, have minimal overhead.

Note that each run yields two distinct sets of MATRIX1 metrics. This is because each execution of run.sh spawns two madevent applications (the "G1" and "G2" subprocesses), for the two MadEvent sampling channels associated to the two Feynman diagrams in this process.

The numbers of MATRIX1 calls indicate the number of phase space points before unweighting. It is puzzling that the ratio between unweighted and weighted events (i.e. the unweighting efficiency) seems to vary: this is equal to approximately 1/8, 1/24 and 1/20 in the three runs. In the gridpack creation step, the unweighting efficiency had been 1.6k/28k, i.e. approximately 1/18.

### Flamegraph (1)

A flamegraph produced with the following script for 100k events
```
  ../flgrAV time ./run.sh 100000 1234
```
is shown below. 

<img src="eemumu/eemumu_madevent_100Kunw-bis.png"  width="1200"/>

The execution took 80s of user CPU time. The sampling rate is set at 1kHz, so 80k frames were collected. 
- Note that the `MATRIX1_` function, where matrix elements (MEs) are computed, only took 4.8k frames, i.e. 4.8s of user CPU time.
- Using the counters to estimate in 2.05M the total number of MATRIX1 calls, this gives an approximate throughput of 2.05M/4.8s i.e. 4.2E5 MEs per second for this Fortran implementation.

## EEMUMU gridpack: time performance for unweighted event generation (2)

Some performance timers have then been added later been added for the Fortran subroutine MATRIX1.

On itscrd70, generating a variable number of unweighted events
```
  time ./run.sh <nunw> 1234
```
gives the following approximate performance (using gcc/gfortran 4.8.5, the default system compiler):

|   nunw |   real |  user |  sys |    MATRIX1 calls |  MATRIX1 times | MATRIX1 throughputs |
|   ---: |   ---: |  ---: | ---: |             ---: |           ---: |                ---: |
|   1000 |   1.9s |  0.9s | 0.2s |       4019, 4019 | 0.013s, 0.013s |  3.04E5/s, 3.04E5/s |
|  10000 |  14.4s | 10.4s | 1.3s |   124019, 124019 | 0.404s, 0.409s |  3.07E5/s, 3.03E5/s |
| 100000 | 103.2s | 83.9s | 9.4s | 1020019, 1020019 | 3.320s, 3.331s |  3.07E5/s, 3.06E5/s |

The timers have introduced a 5-10% overhead to user and sys CPU times, with respect to those in the previous table.

The overhead, however, is much larger for the MATRIX1 calls themselves, which are estimated to take 6.6s overall, instead of the 4.8s previously estimated from the flamegraph. As a consequence, the 3.0E5/s throughputs in this table are probably underestimated.

### Flamegraph (2)

A flamegraph produced with the following script for 100k events
```
  ../flgrAV time ./run.sh 100000 1234
```
is shown below. 

<img src="eemumu/eemumu_madevent_100Kunw-bis2.png"  width="1200"/>

The execution took 86s of user CPU time. The sampling rate is set at 1kHz, so 86k frames were collected. 
- The `MATRIX1_` function, where MEs are computed, took 6.9k frames, i.e. 6.9s of user CPU time. Out of these, 2.1s are an overhead from the counter start/stop functions. The actual ME calculation probably took 4.8s, as previously estimated. Hence, the previous throughput estimate of 4.2E5/s seems correct.
