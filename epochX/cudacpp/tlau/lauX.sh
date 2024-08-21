#!/bin/bash
# Copyright (C) 2020-2023 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jun 2023) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2023) for the MG5aMC CUDACPP plugin.

set -e # fail on error

scrdir=$(cd $(dirname $0); pwd)

function usage()
{
  echo "Usage:   $0 -<backend> <procdir> [-nomakeclean] [-rndoff <integer_offset|x10>] [-togridpack|-fromgridpack]"
  echo "         (supported <backend> values: fortran, cuda, hip, cppnone, cppsse4, cppavx2, cpp512y, cpp512z)"
  echo "         (*special* <backend> value: ALL processes all available backends)"
  echo "Example: $0 -cppavx2 gg_tt.mad"
  exit 1
}

bckend=
proc=
grid=
nomakeclean=
rndoff=0
while [ "$1" != "" ]; do
  if [ "$1" == "-fortran" ] || [ "$1" == "-cuda" ] || [ "$1" == "-hip" ] || [ "$1" == "-cppnone" ] || [ "$1" == "-cppsse4" ] || [ "$1" == "-cppavx2" ] || [ "$1" == "-cpp512y" ] || [ "$1" == "-cpp512z" ] || [ "$1" == "-ALL" ]; then
    if [ "${bckend}" == "" ]; then bckend=${1/-/}; else echo "ERROR! Backend already set"; usage; fi
  elif [ "$1" == "-togridpack" ] && [ "${grid}" == "" ]; then
    grid=${1}
  elif [ "$1" == "-fromgridpack" ] && [ "${grid}" == "" ]; then
    grid=${1}
  elif [ "$1" == "-nomakeclean" ]; then
    nomakeclean=$1
  elif [ "$1" == "-rndoff" ]; then
    rndoff=$2
    shift
  elif [ "${1#-}" != "${1}" ]; then
    echo "ERROR! Invalid option '$1'"
    usage
  elif [ "${proc}" == "" ]; then
    proc=$1
  else
    echo "ERROR! Invalid input '$1': process directory already set to '${proc}'"
    usage
  fi
  shift
done
if [ "${bckend}" == "" ]; then echo "ERROR! No backend was specified"; usage; fi
if [ "${proc}" == "" ]; then echo "ERROR! No process directory was specified"; usage; fi
if [ "${grid}" != "" ] && [ "${rndoff}" != "0" ]; then echo "ERROR! ${grid} and -rndoff are not compatible"; exit 1; fi

if [ "${bckend}" == "ALL" ]; then
  if [ "${grid}" == "-togridpack" ]; then echo "ERROR! ${grid} and -ALL are not compatible"; exit 1; fi # temporary?
  for b in fortran cuda hip cppnone cppsse4 cppavx2 cpp512y cpp512z; do
    $0 -${b} ${nomakeclean} ${proc} -rndoff ${rndoff} ${grid}
    nomakeclean=-nomakeclean # respect user input only on the first test, then keep the builds
  done
  exit 0 # successful termination on each loop (and skip the rest of this file)
fi

gridpackdir=${scrdir}/gridpacks/${proc}
suff=.mad
if [ "${proc}" == "${proc%${suff}}" ]; then echo "ERROR! Process directory does not end in '${suff}'"; usage; fi
proc=${proc%${suff}}
resultsdir=${scrdir}/logs_${proc//_}_${bckend}
if [ "${grid}" == "-togridpack" ]; then
  resultsdir=${gridpackdir}
elif [ "${grid}" == "-fromgridpack" ]; then
  resultsdir=${gridpackdir/gridpacks/fromgridpacks}/${bckend}
  rm -rf ${resultsdir}; mkdir -p ${resultsdir}
fi

if [ "${rndoff}" == "x10" ]; then
  for i in $(seq 0 9); do
    $0 -${bckend} ${nomakeclean} ${proc} ${grid} -rndoff ${i}
    nomakeclean=-nomakeclean # respect user input only on the first test, then keep the builds
  done
  more ${resultsdir}/*txt | \egrep '(Cross-section :)'
  exit 0 # successful termination on each loop (and skip the rest of this file)
elif [ ${rndoff} -lt 0 ]; then
  echo "ERROR! Invalid rndoff=${rndoff}"
  exit 1
fi

outfile=output.txt
if [ "${rndoff}" != "0" ]; then outfile=output${rndoff}.txt; fi

function exit0()
{
  echo ""
  echo "********************************************************************************"
  echo ""
  exit 0
}

if [ "${bckend}" == "cuda" ]; then
  if ! nvidia-smi -L > /dev/null 2>&1; then
    echo "WARNING! No NVidia GPU was found: skip backend ${bckend}"
    exit0
  fi
elif [ "${bckend}" == "hip" ]; then
  if ! rocm-smi -i > /dev/null 2>&1; then
    echo "WARNING! No AMD GPU was found: skip backend ${bckend}"
    exit0
  fi
fi

function getnevt()
{
  if [ "${proc}" == "gg_tt" ]; then
    nevt=10000
  elif [ "${proc}" == "gg_ttg" ]; then
    nevt=10000
  elif [ "${proc}" == "gg_ttgg" ]; then
    nevt=1000
  elif [ "${proc}" == "gg_ttggg" ]; then
    nevt=100
  else
    echo "WARNING! Unknown process ${proc}" > /dev/stderr
    nevt=100
  fi
  echo $nevt
}

# Get the random seed
(( rndseed = 21 + ${rndoff} ))

# Get the number of unweighted events to generate
nevt=$(getnevt)

function showcpugpu()
{
  unames=$(uname -s)
  unamep=$(uname -p)
  if nvidia-smi -L > /dev/null 2>&1; then
    gpuTxt="$(nvidia-smi -L | wc -l)x $(nvidia-smi -L | awk '{print $3,$4}' | sort -u)"
  elif rocm-smi -i > /dev/null 2>&1; then
    gpuTxt="$(rocm-smi --showproductname | grep 'Card series' | awk '{print $5,$6,$7}')"
  else
    gpuTxt=none
  fi
  if [ "${unames}" == "Darwin" ]; then 
    cpuTxt=$(sysctl -h machdep.cpu.brand_string)
    cpuTxt=${cpuTxt/machdep.cpu.brand_string: }
  elif [ "${unamep}" == "ppc64le" ]; then 
    cpuTxt=$(cat /proc/cpuinfo | grep ^machine | awk '{print substr($0,index($0,"Power"))", "}')$(cat /proc/cpuinfo | grep ^cpu | head -1 | awk '{print substr($0,index($0,"POWER"))}')
  else
    cpuTxt=$(cat /proc/cpuinfo | grep '^model name' |& head -1 | awk '{i0=index($0,"Intel"); if (i0==0) i0=index($0,"AMD"); i1=index($0," @"); if (i1>0) {print substr($0,i0,i1-i0)} else {print substr($0,i0)}}')
  fi
  cpuTxt="${cpuTxt} (nproc=$(nproc))"
  echo -e "On $HOSTNAME [CPU: $cpuTxt] [GPU: $gpuTxt]:"
}

# --- OPTION 2: RUN FROM GRIDPACK ---
if [ "${grid}" == "-fromgridpack" ]; then
  echo "Execute $(basename $0) for process ${proc} and backend ${bckend} from gridpack directory $(pwd)"
  if [ ! -d ${gridpackdir} ]; then echo "ERROR! Gridpack directory '${gridpackdir}' does not exist"; usage; fi
  cd ${gridpackdir}
  rm -rf madevent run.sh events.lhe*
  tar -xzf run_01_gridpack.tar.gz
  # Configure gridpack patches
  dir=madevent/bin/internal
  pushd $dir >& /dev/null
  echo "INFO: configure gridpack patches in ${dir}"
  mv madevent_interface.py madevent_interface.py.BKP
  mv gen_ximprove.py gen_ximprove.py.BKP
  mv cluster.py cluster.py.BKP
  \cp ../../../../MG5aMC_patches/madevent_interface.py .
  \cp ../../../../MG5aMC_patches/gen_ximprove.py .
  \cp ../../../../MG5aMC_patches/cluster.py .
  popd >& /dev/null
  # Configure the appropriate backend
  for dir in madevent/SubProcesses/P*; do
    pushd $dir >& /dev/null
    echo "INFO: redefine madevent symlink for backend ${bckend} in ${dir}"
    if [ "${bckend}" == "fortran" ]; then
      exe=$(\ls madevent_fortran)
    else
      exe=$(\ls build.*${bckend/cpp/}*/madevent_*)
    fi
    if [ "$(echo ${exe} | wc -w)" == "0" ]; then
      echo "ERROR! No madevent executable found for backend $bckend"; exit 1
    elif [ "$(echo ${exe} | wc -w)" != "1" ]; then
      echo "ERROR! Too many madevent executables found for backend $bckend:"; echo "$exe"; exit 1
    fi
    \rm -f madevent
    ln -sf $exe madevent
    echo "----> madevent symlink will now point to $exe"
    popd >& /dev/null
  done
  # Run the test for the appropriate backend
  START=$(date +%s)
  echo "START: $(date)" |& tee ${resultsdir}/${outfile}
  showcpugpu |& tee -a ${resultsdir}/${outfile}
  if [ -v CUDACPP_RUNTIME_DISABLEFPE ]; then echo CUDACPP_RUNTIME_DISABLEFPE is set |& tee -a ${resultsdir}/${outfile}; else echo CUDACPP_RUNTIME_DISABLEFPE is not set |& tee -a ${resultsdir}/${outfile}; fi # temporary? (debug FPEs in CMS DY #942)
  if [ -v CUDACPP_RUNTIME_SKIPXBINCHECKS ]; then echo CUDACPP_RUNTIME_SKIPXBINCHECKS is set |& tee -a ${resultsdir}/${outfile}; else echo CUDACPP_RUNTIME_SKIPXBINCHECKS is not set |& tee -a ${resultsdir}/${outfile}; fi
  ls -l madevent/SubProcesses/P*/madevent |& tee -a ${resultsdir}/${outfile}
  ./run.sh ${nevt} ${rndseed} |& tee -a ${resultsdir}/${outfile}
  mv events* ${resultsdir}
  echo "END: $(date)" |& tee -a ${resultsdir}/${outfile}
  END=$(date +%s)
  echo "ELAPSED: $((END-START)) seconds" |& tee -a ${resultsdir}/${outfile}
  exit0
fi
  
# --- OPTION 1: RUN FROM SOURCE CODE (AND OPTIONALLY CREATE A GRIDPACK) ---
cd $(dirname $0)/..
echo "Execute $(basename $0) for process ${proc} and backend ${bckend} in directory $(pwd)"
procdir=$(pwd)/${proc}${suff}
if [ ! -d ${procdir} ]; then echo "ERROR! Process directory '${procdir}' does not exist"; usage; fi
cd ${procdir}

function lauX_makeclean()
{
  if [ "${nomakeclean}" == "" ]; then
    echo "INFO: clean all builds"
    for d in SubProcesses/P*; do cd $d; make cleanall > /dev/null; cd - > /dev/null; break; done
  else
    echo "WARNING! Keep all builds (-nomakeclean option was specified)"
  fi
}

function lauX_cleanup()
{
  rm -f crossx.html index.html
  rm -f SubProcesses/results.dat
  rm -rf Events HTML; mkdir Events HTML; touch Events/.keep HTML/.keep
  for d in SubProcesses/P*; do cd $d; rm -rf gensym input_app.txt symfact.dat G[0-9]* ajob[0-9]*; cd - > /dev/null; done
  if [ "${grid}" == "-togridpack" ]; then
    rm -rf bin/TheChopper-pl
    rm -rf bin/clean4grid
    rm -rf bin/compile
    rm -rf bin/gridrun
    rm -rf bin/internal/gen_ximprove
    rm -rf bin/refine4grid
    rm -rf bin/replace.pl
    rm -rf bin/run.sh
  fi
}

# Clean builds before launch
lauX_makeclean

# Back up config before launch
cp SubProcesses/randinit SubProcesses/randinit.BKP # save the initial file
cp Cards/run_card.dat Cards/run_card.dat.BKP # save the initial file
cp Cards/grid_card.dat Cards/grid_card.dat.BKP # save the initial file
cp Source/run_card.inc Source/run_card.inc.BKP # save the initial file
cp bin/internal/gen_ximprove.py bin_internal_gen_ximprove.py.BKP # save the initial file
cp bin/internal/madevent_interface.py bin_internal_madevent_interface.py.BKP # save the initial file
cp Source/make_opts Source/make_opts.BKP # save the initial file
cp Source/param_card.inc Source/param_card.inc.BKP # save the initial file

# Clean config before launch
# (NB: "just in case" actions below should normally keep the defaults of generated code in the repo?)
# (NB: but some small differences have been observed, e.g. "False     = gridpack" vs "False = gridpack")
if [ "${rndoff}" == "0" ]; then rm -rf ${resultsdir}; mkdir ${resultsdir}; fi
lauX_cleanup
rm -f SubProcesses/ME5_debug
echo "r=21" > SubProcesses/randinit # just in case
sed -i "s/.* = nevents/  10000 = nevents/" Cards/run_card.dat # just in case
sed -i "s/.* = cudacpp_backend/ cpp = cudacpp_backend/" Cards/run_card.dat # just in case
sed -i "s/.* = cudacpp_bldall/ False = cudacpp_bldall/" Cards/run_card.dat # just in case
sed -i "s/.* = gridpack/  False = gridpack/" Cards/run_card.dat # just in case
sed -i "s/      NEVENTS = .*/      NEVENTS = 10000/" Source/run_card.inc # just in case
sed -i "s/8192 1 1/%(event)s         %(maxiter)s           %(miniter)s/" bin/internal/gen_ximprove.py # just in case
sed -i "s/'int', 8192,'Number of points/'int', 1000,'Number of points/" bin/internal/madevent_interface.py # just in case
sed -i "s/'int', 1, 'Number of iterations'/'int', 5, 'Number of iterations'/" bin/internal/madevent_interface.py # just in case

# Set the random seed
echo "r=${rndseed}" > SubProcesses/randinit # just in case a previous test was not cleaned up

# Set the number of events and iterations in the survey step
sed -i "s/'int', 1000,'Number of points/'int', 8192,'Number of points/" bin/internal/madevent_interface.py
sed -i "s/'int', 5, 'Number of iterations'/'int', 1, 'Number of iterations'/" bin/internal/madevent_interface.py

# Set the number of events and iterations in the refine step
sed -i "s/%(event)s         %(maxiter)s           %(miniter)s/8192 1 1/" bin/internal/gen_ximprove.py

# Set the number of unweighted events in run_card.dat
sed -i "s/ 10000 = nevents/ ${nevt} = nevents/" Cards/run_card.dat

# Set the backend in run_card.dat
sed -i "s/ cpp = cudacpp_backend/${bckend} = cudacpp_backend/" Cards/run_card.dat

# Set gridpack mode in run_card.dat
if [ "${grid}" == "-togridpack" ]; then sed -i "s/.* = gridpack/  True = gridpack/" Cards/run_card.dat; fi

# Configure bldall in run_card.dat
sed -i "s/.* = cudacpp_bldall/ True = cudacpp_bldall/" Cards/run_card.dat

# Launch (generate_events)
# (BUG #683: generate_events does not return an error code even if it fails)
###set -x # verbose
START=$(date +%s)
echo "START: $(date)" |& tee ${resultsdir}/${outfile}
showcpugpu |& tee -a ${resultsdir}/${outfile}
if [ -v CUDACPP_RUNTIME_DISABLEFPE ]; then echo CUDACPP_RUNTIME_DISABLEFPE is set |& tee -a ${resultsdir}/${outfile}; else echo CUDACPP_RUNTIME_DISABLEFPE is not set |& tee -a ${resultsdir}/${outfile}; fi # temporary? (debug FPEs in CMS DY #942)
MG5AMC_CARD_PATH=$(pwd)/Cards time ./bin/generate_events -f |& tee -a ${resultsdir}/${outfile}
echo "END: $(date)" |& tee -a ${resultsdir}/${outfile}
END=$(date +%s)
echo "ELAPSED: $((END-START)) seconds" |& tee -a ${resultsdir}/${outfile}
###set +x # not verbose

# Copy output gridpack to tlau/gridpacks directory
if [ "${grid}" == "-togridpack" ]; then
  mv run_01_gridpack.tar.gz ${gridpackdir}/run_01_gridpack.tar.gz
  echo "Gridpack created: ${gridpackdir}/run_01_gridpack.tar.gz"
fi

# Process and keep results (only for the default rndoff)
if [ "${grid}" != "-togridpack" ] && [ "${rndoff}" == "0" ]; then
  \rm HTML/results.pkl
  mv Events ${resultsdir}; mv HTML ${resultsdir}
  gunzip ${resultsdir}/Events/run_01/unweighted_events.lhe.gz
  # FIXME! No need to keep events in git, there is no lhe file comparison yet anyway (20-DEC-2023)
  \rm ${resultsdir}/Events/run_01/unweighted_events.lhe
  \rm ${resultsdir}/Events/run_01/run_01_tag_1_banner.txt
  touch ${resultsdir}/Events/run_01/.keep
fi

# Clean config after launch
lauX_cleanup
mv SubProcesses/randinit.BKP SubProcesses/randinit # restore the initial file
mv Cards/run_card.dat.BKP Cards/run_card.dat # restore the initial file
mv Cards/grid_card.dat.BKP Cards/grid_card.dat # restore the initial file
mv Source/run_card.inc.BKP Source/run_card.inc # restore the initial file
mv bin_internal_gen_ximprove.py.BKP bin/internal/gen_ximprove.py # restore the initial file
mv bin_internal_madevent_interface.py.BKP bin/internal/madevent_interface.py # restore the initial file
mv Source/make_opts.BKP Source/make_opts # restore the initial file
mv Source/param_card.inc.BKP Source/param_card.inc # restore the initial file

# Add an 80-character separator and exit
exit0
