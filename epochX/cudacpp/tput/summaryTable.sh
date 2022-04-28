#!/bin/sh

table=
if [ "$1" == "-ALL" ] && [ "$2" == "" ]; then
  $0
  $0 -hrdcod
  $0 -juwels
  $0 -alpaka
  $0 -macm1
  $0 -alphas
  exit 0
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
fi
if [ "$1" != "" ]; then echo "Usage: $0 [-ALL|-hrdcod|-juwels|-alpaka|-macm1|-alphas]"; exit 1; fi

unames=$(uname -s)
if [ "${unames}" == "Darwin" ]; then
  echo "ERROR! This script is not supported on MacOS (issues in awk): please execute it on Linux"
  exit 1
fi

cd $(dirname $0)/..
echo PWD=$(pwd)

# Output file
if [ "$table" == "" ]; then
  out=tput/summaryTable.txt
else
  out=tput/summaryTable_${table}.txt
fi
\rm -f $out
touch $out

#----------------------------------------------------------------------------
# A few reference points 
# 1. Code (21 Apr 2022, c0c27684): just before the alphas PR #434
#    Logs (28 Apr 2022, 88fe36d1): cuda116/gcc102 (56 logs from allTees.sh)
# 2. Code (28 Apr 2022, bae5c248): about to be merged in the alphas PR #434
#    Logs (28 Apr 2022, bae5c248): cuda116/gcc102 (56 logs from allTees.sh)
#----------------------------------------------------------------------------

# Select revisions of cudacpp and alpaka logs
crevs=""
arevs=""
if [ "$table" == "" ]; then
  crevs="$crevs 09e482e"  # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  crevs="$crevs 88dc717"  # cuda116/icx2022 (03 Mar 2022) ICX TEST eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
elif [ "$table" == "hrdcod" ]; then
  crevs="$crevs 09e482e"  # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  crevs="$crevs 88dc717"  # cuda116/icx2022 (03 Mar 2022) ICX TEST eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
elif [ "$table" == "juwels" ]; then
  crevs="$crevs 09e482e"  # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  crevs="$crevs 65730b2"  # cuda115/gcc112  (18 Feb 2022) JUWELSCL eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
  crevs="$crevs df441ad"  # cuda115/gcc112  (18 Feb 2022) JUWELSBO eemumu/ggtt/ggttgg x f/d x inl0/inl1 + ggttg/ggttggg x f/d
elif [ "$table" == "alpaka" ]; then
  crevs="$crevs 09e482e"  # cuda116/gcc102  (03 Mar 2022) BASELINE eemumu/ggtt/ggttgg x f/d x hrd0/hrd1 x inl0/inl1 + ggttg/ggttggg x f/d x hrd0/hrd1 x inl0
  arevs="$arevs b93a2f3"  # cuda116/gcc102  (06 Mar 2022) GOLDEPX4 eemumu/ggtt/ggttg/ggttgg/ggttggg x d x inl0
elif [ "$table" == "macm1" ]; then
  crevs="$crevs ea28661"  # NOGPU/clang120  (10 Apr 2022) MACM1ARM eemumu/ggtt/ggttg/ggttgg/ggttggg x f/d x hrd0 x inl0
elif [ "$table" == "alphas" ]; then
  crevs="$crevs 88fe36d1" # cuda116/gcc102  (28 Apr 2022) PRE-AS   56 logs allTees.sh
  crevs="$crevs bae5c248" # cuda116/gcc102  (28 Apr 2022) POST-AS  56 logs allTees.sh
fi

# Select processes
procs="eemumu ggtt ggttg ggttgg ggttggg"

# Select fptype, helinl, hrdcod
if [ "$table" == "" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
elif [ "$table" == "hrdcod" ]; then
  fpts="d"
  inls="inl0"
  hrds="hrd0 hrd1"
elif [ "$table" == "juwels" ]; then
  fpts="d f"
  inls="inl0 inl1"
  hrds="hrd0"
elif [ "$table" == "alpaka" ]; then
  fpts="d"
  inls="inl0"
  hrds="hrd0"
elif [ "$table" == "macm1" ]; then
  fpts="d f"
  inls="inl0"
  hrds="hrd0"
elif [ "$table" == "alphas" ]; then
  fpts="d f"
  inls="inl0" # no need to add inl1 (even if data exists)
  hrds="hrd0 hrd1"
fi

# Select tag list
if [ "$table" == "alpaka" ]; then
  taglist="CUD/none ALP/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z"
elif [ "$table" == "macm1" ]; then
  taglist="CPP/none CPP/sse4"
else
  taglist="CUD/none CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z"
fi

# Kernel function
function oneTable()
{
  files=""
  for proc in $procs; do
    file=../${bckend}/tput/logs_${proc}_${suff}/log_${proc}_${suff}_${fpt}_${inl}_${hrd}.txt
    if [ -f $file ]; then files="$files $file"; fi
  done
  ###echo "*** FILES $files ***" >> $out
  if [ "$files" == "" ]; then continue; fi
  git checkout $rev $files >& /dev/null
  node=$(cat $files | grep ^On | sort -u)
  if [ "$nodelast" != "$node" ]; then echo -e "$node\n" >> $out; nodelast=$node; fi
  ###cat $files | awk '/^runExe.*check.*/{print $0};/^Process/{print $0};/Workflow/{print $0};/MECalcOnly/{print $0}'; continue
  cat $files | awk -vtaglist="$taglist" -vrev=$rev -vcomplast=none -vinllast=none -vhrdlast=none -vfptlast=none '\
    /^runExe .*check.*/{split($0,a,"check.exe"); last=substr(a[1],length(a[1])); if (last=="g"){tag="CUD"} else if(last=="p"){tag="ALP"} else{tag="CPP"}; split($0,a,"build."); split(a[2],b,"_"); tag=tag"/"b[1]};\
    /^runExe .*check.*/{split($0,a," -p "); split(a[2],b); grid=b[1]"/"b[2]"/"b[3]};\
    /^runExe .*check.*/{proc=$0; gsub("/SubProcesses.*","",proc); gsub(".*/","",proc); gsub(".auto","",proc); grid_proc_tag[proc,tag]=grid};\
    ###/^Process/{split($3,a,"_"); proc=a[3]"_"a[4]; grid_proc_tag[proc,tag]=grid};\
    /^Process(.)/{split($0,a,"["); comp="["a[2]; if ( complast == "none" ){print comp; complast=comp}};\
    /^Process/{split($0,a,"]"); split(a[2],b,"="); inl=b[2]; if(inl!=inllast){printf "HELINL="inl; inllast=inl}}\
    /^Process/{split($0,a,"]"); split(a[3],b,"="); hrd=b[2]; if(hrd!=hrdlast){if(hrd==""){hrd=0}; print " HRDCOD="hrd; hrdlast=hrd}}\
    /^FP precision/{fpt=$4; /*if ( fpt != fptlast ){print "FPTYPE="fpt; fptlast=fpt}*/}\
    ###/Workflow/{split($4,a,":"); tag=a[1]; split($4,a,"+"); split(a[4],b,"/"); tag=tag"/"b[2]};\
    /.*check.exe: Aborted/{tput_proc_tag[proc,tag]="(FAILED)"};\
    /MECalcOnly/{tput=sprintf("%.2e", $5); tput_proc_tag[proc,tag]=tput};\
    END{ntag=split(taglist,tags);\
        ###nproc=split("EPEM_MUPMUM GG_TTX GG_TTXG GG_TTXGG GG_TTXGGG",procs);\
        ###procs_txt["EPEM_MUPMUM"]="eemumu";\
        ###procs_txt["GG_TTX"]="ggtt";\
        ###procs_txt["GG_TTXG"]="ggttg";\
        ###procs_txt["GG_TTXGG"]="ggttgg";\
        ###procs_txt["GG_TTXGGG"]="ggttggg";\
        nproc=split("ee_mumu gg_tt gg_ttg gg_ttgg gg_ttggg",procs);\
        procs_txt["ee_mumu"]="eemumu";\
        procs_txt["gg_tt"]="ggtt";\
        procs_txt["gg_ttg"]="ggttg";\
        procs_txt["gg_ttgg"]="ggttgg";\
        procs_txt["gg_ttggg"]="ggttggg";\
        printf "%8s", ""; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; printf "%14s", procs_txt[proc]}; printf "\n";\
        gridslast="";\
        for(itag=1;itag<=ntag;itag++)\
        {tag=tags[itag];\
         gridsok=0; grids=""; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; grid=grid_proc_tag[proc,tag]; if(grid==""){grid="------"} else {gridsok=1}; grids=grids""sprintf("%14s","["grid"]")}; grids=grids"\n";\
         if(grids!=gridslast && gridsok==1){printf "%-8s%s", "", grids}; gridslast=grids;\
         printf "%-8s", tag; for(iproc=1;iproc<=nproc;iproc++){proc=procs[iproc]; tput=tput_proc_tag[proc,tag]; if(tput==""){tput="--------"}; printf "%14s", tput}; printf "\n";\
         }}' >> $out
  echo "" >> $out
}
          
# Iterate through log files
for fpt in $fpts; do
  echo -e "*** FPTYPE=$fpt ******************************************************************\n" >> $out
  for bckend in cudacpp alpaka; do
    if [ "$bckend" == "alpaka" ]; then
      revs="$arevs"
      suff=auto
    else
      revs="$crevs"
      suff=manu
    fi
    if [ "$revs" == "" ]; then continue; fi
    ### New sorting (alphas)
    if [ "$table" == "alphas" ]; then
      for inl in $inls; do
        for hrd in $hrds; do
          echo -e "-------------------------------------------------------------------------------\n" >> $out
          for rev in $revs; do
            echo -e "+++ $bckend REVISION $rev +++" >> $out
            nodelast=
            oneTable
          done
        done
      done
    else
      ### Old sorting (all but alphas)
      for rev in $revs; do
        echo -e "+++ $bckend REVISION $rev +++" >> $out
        nodelast=
        for inl in $inls; do
          for hrd in $hrds; do
            oneTable
          done
        done
      done
    fi
  done
done
git checkout HEAD tput/logs* >& /dev/null
cat $out
