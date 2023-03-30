#!/bin/sh

# NB: some PRs include many log modifications but eventually go back to the initial logs (overall noop)
# NB: use 'git log --full-history' to show all intermediate commits (which are pruned by a simple 'git log')

table=
if [ "$1" == "-ALL" ] && [ "$2" == "" ]; then
  set -e
  $0 -latest
  $0 -bridge
  $0 -hrdcod
  $0 -juwels
  $0 -alpaka
  $0 -macm1
  $0 -alphas
  $0 -3xcomp
  $0 -ol23silv
  ###$0 -ol23gold
  exit 0
elif [ "$1" == "-latest" ]; then
  table="latest"; shift
elif [ "$1" == "-bridge" ]; then
  table="bridge"; shift
elif [ "$1" == "-hrdcod" ]; then
  table="hrdcod"; shift
elif [ "$1" == "-juwels" ]; then
  table="juwels"; shift
elif [ "$1" == "-alpaka" ]; then
  table="alpaka"; shift
elif [ "$1" == "-macm1" ]; then
  table="macm1"; shift
elif [ "$1" == "-alphas" ]; then
  table="alphas"; shift
elif [ "$1" == "-3xcomp" ]; then
  table="3xcomp"; shift
elif [ "$1" == "-ol23silv" ]; then
  table="ol23silv"; shift
elif [ "$1" == "-ol23gold" ]; then
  table="ol23gold"; shift
else
  echo "Usage: $0 <table [-ALL|-latest|-bridge|-hrdcod|-juwels|-alpaka|-macm1|-alphas|-3xcomp|-ol23silv|ol23gold]>"; exit 1
fi

unames=$(uname -s)
if [ "${unames}" == "Darwin" ]; then
  echo "ERROR! This script is not supported on MacOS (issues in awk): please execute it on Linux"
  exit 1
fi

cd $(dirname $0)/..
echo PWD=$(pwd)

# Output file
out=tput/summaryTable_${table}.txt
\rm -f $out
touch $out

#-------------------------------------------------------------------------------------------------------
# A few reference points 
#-------------------------------------------------------------------------------------------------------
# 1. Code (21 Apr 2022, c0c27684): just before the alphas PR #434
#    Logs (28 Apr 2022, 88fe36d1): cuda116/gcc102 (56 logs from allTees.sh)
# 2. Code (28 Apr 2022, bae5c248): about to be merged in the alphas PR #434
#    Logs (28 Apr 2022, bae5c248): cuda116/gcc102 (56 logs from allTees.sh)
#-------------------------------------------------------------------------------------------------------
# 1. Logs (20 Jun 2022, 2d3e789): cuda116/gcc102 (60 logs from allTees.sh) <= previous cuda116/gcc102
# 2. Logs (26 Aug 2022, 21c4cb8): cuda116/gcc102 (60 logs from allTees.sh) <= LATEST cuda116/gcc102
# 3. Logs (26 Aug 2022, d250d2d): cuda117/gcc112 (60 logs from allTees.sh) <= LATEST cuda117/gcc112
#----------------------------------------------------------------------------

# Select revisions of cudacpp and alpaka logs
crevs="" # cudacpp .sa
arevs="" # alpaka
mrevs="" # cudacpp .mad
if [ "$table" == "latest" ]; then
  mrevs="$mrevs d250d2d" # cuda117/gcc112  (26 Aug 2022) GCC112   60 logs allTees.sh
elif [ "$table" == "bridge" ]; then
  mrevs="$mrevs d250d2d" # cuda117/gcc112  (26 Aug 2022) GCC112   60 logs allTees.sh
elif [ "$table" == "hrdcod" ]; then
  mrevs="$mrevs d250d2d" # cuda117/gcc112  (26 Aug 2022) GCC112   60 logs allTees.sh
  mrevs="$mrevs f3ee68c" # cuda117/icx2022 (27 Aug 2022) ICX2022  60 logs allTees.sh
elif [ "$table" == "juwels" ]; then
  crevs="$crevs 09e482e" # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  crevs="$crevs 65730b2" # cuda115/gcc112  (18 Feb 2022) JUWELSCL eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
  crevs="$crevs df441ad" # cuda115/gcc112  (18 Feb 2022) JUWELSBO eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
  ###mrevs="$mrevs 2d3e789" # cuda116/gcc102  (20 Jun 2022) BASELINE eemumu/ggtt* x f/d x hrd0 x inl0 x default/bridge
  ###mrevs="$mrevs fd4ffe1" # cuda115/gcc112  (24 Jun 2022) JUWELSCL eemumu/ggtt* x f/d x hrd0 x inl0 x default/bridge
elif [ "$table" == "alpaka" ]; then
  crevs="$crevs 09e482e" # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  arevs="$arevs f5a44ba" # cuda116/gcc102  (06 Mar 2022) GOLDEPX4 eemumu/ggtt/ggttg/ggttgg/ggttggg x d x inl0 (golden tag epochx4)
elif [ "$table" == "macm1" ]; then
  crevs="$crevs ea28661" # NOGPU/clang120  (10 Apr 2022) MACM1ARM eemumu/ggtt/ggttg/ggttgg/ggttggg x f/d x hrd0 x inl0
elif [ "$table" == "alphas" ]; then
  crevs="$crevs 88fe36d" # cuda116/gcc102  (28 Apr 2022) PRE-AS   56 logs allTees.sh
  crevs="$crevs bae5c24" # cuda116/gcc102  (28 Apr 2022) POST-AS  56 logs allTees.sh
elif [ "$table" == "3xcomp" ]; then
  mrevs="$mrevs d250d2d" # cuda117/gcc112  (26 Aug 2022) GCC112   60 logs allTees.sh
  mrevs="$mrevs 2c77b32" # cuda117/clang13 (27 Aug 2022) CLANG13  60 logs allTees.sh
  mrevs="$mrevs f3ee68c" # cuda117/icx2022 (27 Aug 2022) ICX2022  60 logs allTees.sh
elif [ "$table" == "ol23silv" ]; then
  ###mrevs="$mrevs 49e4381" # cuda120/gcc112  (23 Feb 2023 itscrd90)
  mrevs="$mrevs 2a6ddd0" # cuda120/gcc121  (24 Feb 2023 itscrd90)
  mrevs="$mrevs aed7c84" # cuda120/clang14 (23 Feb 2023 itscrd90)
  mrevs="$mrevs 32e098f" # cuda120/icx2023 (23 Feb 2023 itscrd90)
elif [ "$table" == "ol23gold" ]; then
  mrevs="$mrevs e67558d" # nocuda/gcc121   (10 Mar 2023 bmkGold6130)
  mrevs="$mrevs aa4a5c4" # nocuda/clang14  (09 Mar 2023 bmkGold6130)
  mrevs="$mrevs 99c1bd6" # nocuda/icx2023  (10 Mar 2023 bmkGold6130)
else
  echo "ERROR! Unknown table '$table'"; exit 1
fi

# Select processes
procs="eemumu ggtt ggttg ggttgg ggttggg"

# Select fptype, helinl, hrdcod
if [ "$table" == "latest" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
  brds="nobr"
elif [ "$table" == "bridge" ]; then
  fpts="d f"
  inls="inl0"
  hrds="hrd0"
  brds="nobr brdg"
elif [ "$table" == "hrdcod" ]; then
  fpts="d"
  inls="inl0"
  hrds="hrd0 hrd1"
  brds="nobr"
elif [ "$table" == "juwels" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
  ###brds="nobr"
  brds="nobr brdg"
elif [ "$table" == "alpaka" ]; then
  fpts="d"
  inls="inl0"
  hrds="hrd0"
  brds="nobr"
elif [ "$table" == "macm1" ]; then
  fpts="d f"
  inls="inl0"
  hrds="hrd0"
  brds="nobr"
elif [ "$table" == "alphas" ]; then
  fpts="d f"
  inls="inl0" # no need to add inl1 (even if data exists)
  hrds="hrd0 hrd1"
  brds="nobr"
elif [ "$table" == "3xcomp" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0 hrd1"
  brds="nobr"
elif [ "$table" == "ol23silv" ] || [ "$table" == "ol23gold" ]; then
  ###fpts="d f m"
  fpts="d f"
  inls="inl0 inl1"
  ###hrds="hrd0 hrd1"
  hrds="hrd0"
  brds="nobr"
else
  echo "ERROR! Unknown table '$table'"; exit 1
fi

# Select tag list
if [ "$table" == "alpaka" ]; then
  taglist="CUD/none ALP/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z"
elif [ "$table" == "macm1" ]; then
  taglist="CPP/none CPP/sse4"
elif [ "$table" == "bridge" ]; then
  taglist="CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z CUD/none" # revert order, CUDA last
else
  taglist="CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z"
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
    ###/^runExe .*check.*/{print $0}
    /^runExe .*check.*/{split($0,a,"check.exe"); last=substr(a[1],length(a[1])); if (last=="g"){tag="CUD"} else if(last=="p"){tag="ALP"} else{tag="CPP"}; split($0,a,"build."); split(a[2],b,"_"); tag=tag"/"b[1]};
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
  for bckend in cudacpp.mad cudacpp.sa alpaka; do
    if [ "$bckend" == "alpaka" ]; then
      revs="$arevs"
      suff=auto
    elif [ "$bckend" == "cudacpp.sa" ]; then
      revs="$crevs"
      suff=manu
    else
      revs="$mrevs"
      suff=mad
    fi
    if [ "$revs" == "" ]; then continue; fi
    ### DIFFERENT SORTINGS
    if [ "$table" == "3xcomp" ] || [ "$table" == "ol23silv" ] || [ "$table" == "ol23gold" ]; then
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
    ### New sorting (alphas)
    elif [ "$table" == "alphas" ]; then
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
    else
      ### Old sorting (all but alphas, 3xcomp, ol23silv, ol23gold)
      for rev in $revs; do
        echo -e "+++ $bckend REVISION $rev (commit date: $(git log $rev --pretty=format:'%ci' --abbrev-commit -n1)) +++" >> $out
        nodelast=
        for inl in $inls; do
          for hrd in $hrds; do
            for brd in $brds; do
              oneTable
            done
          done
        done
      done
    fi
  done
done
###echo TEST >> ${out}
cp ${0} ${0}.NEW
cp ${out} ${out}.NEW
git reset --hard HEAD >& /dev/null
mv ${0}.NEW ${0}
mv ${out}.NEW ${out}
cat $out
