#!/bin/sh
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

# NB: some PRs include many log modifications but eventually go back to the initial logs (overall noop)
# NB: use 'git log --full-history' to show all intermediate commits (which are pruned by a simple 'git log')

echo "==============================================================================="
echo "Executing $0 $*"
echo "==============================================================================="
echo

quiet=
if [ "$1" == "-quiet" ]; then
  quiet=$1; shift
fi

table=
if [ "$1" == "-ALL" ] && [ "$2" == "" ]; then
  set -e
  $0 ${quiet} -rd90_chep24
  $0 ${quiet} -gold_chep24
  $0 ${quiet} -lumi_chep24
  exit 0
elif [ "$1" == "-rd90_chep24" ]; then
  table="rd90_chep24"; shift
elif [ "$1" == "-gold_chep24" ]; then
  table="gold_chep24"; shift
elif [ "$1" == "-lumi_chep24" ]; then
  table="lumi_chep24"; shift
fi

if [ "$1" != "" ] || [ "$table" == '' ]; then
  echo "Usage: $0 [-quiet] <table [-ALL|-rd90_chep24|gold_chep24|lumi_chep24]>"; exit 1
fi

unames=$(uname -s)
if [ "${unames}" == "Darwin" ]; then
  echo "ERROR! This script is not supported on MacOS (issues in awk): please execute it on Linux"
  exit 1
fi

cd $(dirname $0)/..
scr=tput/$(basename ${0})
###echo PWD=$(pwd)

# Output file
out=tput/summaryTable_${table}.txt
\rm -f $out
touch $out

#-------------------------------------------------------------------------------------------------------
# A few reference points 
#-------------------------------------------------------------------------------------------------------
# 1. Logs (03 Oct 2024, cd8e872): cuda120/gcc113 (102 logs from allTees.sh) <= v1.00.00     itscrd90
# 2. Logs (03 Oct 2024, a3d64bd): -------/gcc114 ( 96 logs from allTees.sh) <= v1.00.00     itgold91
# 3. Logs (03 Oct 2024, 07c2a53): roc60/clang170 ( 96 logs from allTees.sh) <= v1.00.00+amd lumi
#----------------------------------------------------------------------------

# Default processes
procs="eemumu ggtt ggttg ggttgg ggttggg"

# Default fptype, helinl, hrdcod, bridge
fpts="d m f" # default
inls="inl0" # default
hrds="hrd0" # default
brds="nobr" # default

# Select revisions of cudacpp logs and modify defaults
if [ "$table" == "rd90_chep24" ]; then
  mrevs="$mrevs cd8e872" # cuda120/gcc113  (03 Oct 2024 itscrd90) v1.00.00
  ###inls="inl0 inl1"
  ###hrds="hrd0 hrd1"
elif [ "$table" == "gold_chep24" ]; then
  mrevs="$mrevs a3d64bd" # -------/gcc114  (03 Oct 2024 gold91)   v1.00.00e
elif [ "$table" == "lumi_chep24" ]; then
  mrevs="$mrevs 07c2a53" # roc60/clang170  (03 Oct 2024 lumi)     v1.00.00+amd
else
  echo "ERROR! Unknown table '$table' (while choosing revisions and modifying defaults)"; exit 1
fi

# Select tag list
if [ "$table" == "macm1" ]; then
  taglist="CPP/none CPP/sse4"
else
  taglist="CUD/none HIP/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z"
fi

# Kernel function
function oneTable()
{
  files=""
  if [ "$brd" == "brdg" ]; then brdsuf="_bridge"; else brdsuf=""; fi
  if [ "${bckend/.*}" == "$(basename $(pwd))" ]; then bckdir=""; else bckdir=../${bckend/.*}/; fi
  for proc in $procs; do
    file=${bckdir}tput/logs_${proc}_${suff}/log_${proc}_${suff}_${fpt}_${inl}_${hrd}${brdsuf}.txt
    if [ "$file" == "$(git ls-tree --name-only $rev $file)" ]; then files="$files $file"; fi
  done
  ###echo "*** FILES $files ***" >> $out
  if [ "$files" == "" ]; then return; fi
  ###echo "git checkout $rev $files"
  git checkout $rev $files >& /dev/null
  if [ "$?" != "0" ]; then echo "ERROR! 'git checkout $rev' failed!"; exit 1; fi
  if [ "$files" == "" ]; then return; fi
  node=$(cat $files | grep ^On | sort -u)
  if [ "$nodelast" != "$node" ]; then echo -e "$node\n" >> $out; nodelast=$node; fi
  ###cat $files | awk '/^runExe.*check.*/{print $0};/^Process/{print $0};/Workflow/{print $0};/MECalcOnly/{print $0}'; continue
  cat $files | awk -vtaglist="$taglist" -vrev=$rev -vcomplast=none -vinllast=none -vhrdlast=none -vfptlast=none -vbrdlast=none '
    /^runExe .*check_cpp.*/{tag="CPP"};
    /^runExe .*check_cpp.*/{split($0,a,"build."); split(a[2],b,"_"); tag=tag"/"b[1]};
    /^runExe .*check_cuda.*/{tag="CUD/none"};
    /^runExe .*check_hip.*/{tag="HIP/none"};
    ###/^runExe .*check.*/{print $0, tag};
    /^runExe .*check.*/{split($0,a," -p "); split(a[2],b); grid=b[1]"/"b[2]"/"b[3]};
    /^runExe .*check.*/{proc=$0; gsub("/SubProcesses.*","",proc); gsub(".*/","",proc); gsub(".auto","",proc); gsub(".mad","",proc); grid_proc_tag[proc,tag]=grid};
    ###/^Process/{split($3,a,"_"); proc=a[3]"_"a[4]; grid_proc_tag[proc,tag]=grid};
    /^Process(.)/{split($0,a,"["); comp="["a[2]; if ( complast == "none" ){print comp; complast=comp}};
    /^Process/{split($0,a,"]"); split(a[2],b,"="); inl=b[2]; if(inl!=inllast){printf "HELINL="inl; inllast=inl}}
    /^Process/{split($0,a,"]"); split(a[3],b,"="); hrd=b[2]; if(hrd!=hrdlast){if(hrd==""){hrd=0}; printf " HRDCOD="hrd; hrdlast=hrd}}
    /Workflow .*+BRD/{brd="yes"; if(brd!=brdlast){print " BRIDGE="brd; brdlast=brd}};
    /Workflow .*+MES/{brd="no"; if(brd!=brdlast){print " BRIDGE="brd; brdlast=brd}};
    ###/Workflow/{split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};
    /^FP precision/{fpt=$4; /*if ( fpt != fptlast ){print "FPTYPE="fpt; fptlast=fpt}*/}
    ###/MECalcOnly/{print proc, tag, $0}
    /MECalcOnly/{tput=sprintf("%.2e", $5); tput_proc_tag[proc,tag]=tput};
    /.*check.exe: Aborted/{tput_proc_tag[proc,tag]="(FAILED)"};
    END{ntag=split(taglist,tags);
        ###nproc=split("EPEM_MUPMUM GG_TTX GG_TTXG GG_TTXGG GG_TTXGGG",procs);
        ###procs_txt["EPEM_MUPMUM"]="eemumu";
        ###procs_txt["GG_TTX"]="ggtt";
        ###procs_txt["GG_TTXG"]="ggttg";
        ###procs_txt["GG_TTXGG"]="ggttgg";
        ###procs_txt["GG_TTXGGG"]="ggttggg";
        nproc=split("ee_mumu gg_tt gg_ttg gg_ttgg gg_ttggg",procs);
        procs_txt["ee_mumu"]="eemumu";
        procs_txt["gg_tt"]="ggtt";
        procs_txt["gg_ttg"]="ggttg";
        procs_txt["gg_ttgg"]="ggttgg";
        procs_txt["gg_ttggg"]="ggttggg";
        printf "%8s", ""; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; printf "%14s", procs_txt[proc]}; printf "\n";
        gridslast="";
        for(itag=1;itag<=ntag;itag++)
        {tag=tags[itag];
         gridsok=0; grids=""; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; grid=grid_proc_tag[proc,tag]; if(grid==""){grid="------"} else {gridsok=1}; grids=grids""sprintf("%14s","["grid"]")}; grids=grids"\n";
         if(grids!=gridslast && gridsok==1){printf "%-8s%s", "", grids}; gridslast=grids;
         printf "%-8s", tag; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; tput=tput_proc_tag[proc,tag]; if(tput==""){tput="--------"}; printf "%14s", tput}; printf "\n";
         }}' >> $out
  echo "" >> $out
}
          
# Iterate through log files
for fpt in $fpts; do
  echo -e "*** FPTYPE=$fpt ******************************************************************\n" >> $out
  for bckend in cudacpp.mad; do
    if [ "$bckend" == "cudacpp.mad" ]; then
      revs="$mrevs"
      suff=mad
    else
      echo "ERROR! Unknown backend '$bckend' (while setting revs and suff)"; exit 1
    fi
    if [ "$revs" == "" ]; then continue; fi
    if [ "${table/chep24}" != "${table}" ]; then
      ### New sorting (3xcomp)
      for inl in $inls; do
        for hrd in $hrds; do
          for brd in $brds; do
            echo -e "-------------------------------------------------------------------------------\n" >> $out
            for rev in $revs; do
              echo -e "+++ $bckend REVISION $rev (commit date: $(git log $rev --pretty=format:'%ci' --abbrev-commit -n1)) +++" >> $out
              nodelast=
              oneTable
            done
          done
        done
      done
    fi
  done
done
###echo TEST >> ${out}
cp ${scr} ${scr}.NEW # allows using a summaryTable.sh not committed to git
cp ${out} ${out}.NEW
git reset --hard HEAD >& /dev/null
mv ${scr}.NEW ${scr} # allows using a summaryTable.sh not committed to git
mv ${out}.NEW ${out}
if [ "$quiet" == "" ]; then cat $out; fi
