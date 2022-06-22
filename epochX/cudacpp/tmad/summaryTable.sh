#!/bin/sh

# Kernel function
function oneTable()
{
  taglist="FORTRAN CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z CUDA/32 CUDA/max"
  parlist="(1) (2-none) (2-sse4) (2-avx2) (2-512y) (2-512z) (3) (3bis)"
  faclist="1 10 100"
  echo "" > $out
  for proc in $procs; do
    file=tmad/logs_${proc}_${suff}/log_${proc}_${suff}_${fpt}_${inl}_${hrd}.txt
    if [ ! -f $file ]; then continue; fi
    ###echo "*** FILE $file ***"
    cat $file | awk -vproc=$proc -vtaglist="$taglist" -vparlist="$parlist" -vfaclist="$faclist" '
      BEGIN{ntag=split(taglist,tags); npar=split(parlist,pars);
            if(ntag!=npar){print "ERROR! ntag!=npar", ntag, npar; exit(1)};
            for(i=1;i<=ntag;i++){tag1[pars[i]]=tags[i]}}
      BEGIN{nfac=split(faclist,facs)}
      BEGIN{lsep=sprintf("%0128d",0); gsub("0","-",lsep)}
      ###/create events.lhe/{print $0}
      /create events.lhe/{par=$2; tag=tag1[par]} # current tag
      /create events.lhe/{fac=substr($5,2)} # current fac
      /\[COUNTERS\]/{if($3=="TOTAL") typ=1; else if($3=="Overhead") typ=2; else if($3=="MEs") typ=3;
                     else{print "ERROR! Unknown type $3"; exit 1};
                     if($4==":") sec=$5; else sec=$8; sec=substr(sec,1,length(sec)-1); # current sec
                     sec3[tag,fac,typ]=sec}
      /\[COUNTERS\]/{if(tag!="" && $3=="MEs")
                     {nevt=$10; ###print tag, nevt;
                      if(tag=="FORTRAN") nevt1[fac]=nevt; 
                      else if(tag=="CUDA/max") nevt1b[fac]=nevt;
                      else if(nevt1[fac]!=nevt){print "ERROR! nevt mismatch", nevt1[fac], nevt; exit 1};
                      tputm1[tag]=tolower($15)}}
      /CHECK/{if(tag!="" && $8=="--bridge"){sabp1[tag]=$5"*"$6"*"$7; sabg1[tag]=$5*$6*$7; gcheck="bridge"}}
      /CHECK/{if(tag!="" && $8!="--bridge"){sap1[tag]=$5"*"$6"*"$7; sag1[tag]=$5*$6*$7; gcheck="nobridge"}}
      /EvtsPerSec/{if(tag!="" && gcheck=="bridge"){tputb1[tag]=$5}}
      /EvtsPerSec/{if(tag!="" && gcheck!="bridge"){tput1[tag]=$5}}
      END{print lsep;
          printf "| %-10s | %-78s | %-30s |\n", proc, "[sec] Total = Overhead (FORTRAN) + MEs (FORTRAN, CPP or CUDA)", "[MEs/sec]"; 
          print lsep;
          for (itag=1; itag<=ntag; itag++)
          {tag=tags[itag]; 
           if(tag=="FORTRAN"){printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "", "mad", "mad", "mad", "sa/brdg", "sa/full";
                              printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "nevt/grid", "32", "32", "32", sabg1["CUDA/32"], sag1["CUDA/32"];
                              printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "nevt total", "x"facs[1]" ["nevt1[facs[1]]"]", "x"facs[2]" ["nevt1[facs[2]]"]", "x"facs[3]" ["nevt1[facs[3]]"]",
                              sabp1["CUDA/32"], sap1["CUDA/32"];
                              print lsep}
           else if(tag=="CUDA/max"){
                              print lsep;
                              printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "", "mad", "mad", "mad", "sa/brdg", "sa/full";
                              printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "nevt/grid", "8192", "8192", "8192", sabg1[tag], sag1[tag];
                              printf "| %-10s | %24s | %24s | %35s | %8s | %8s |\n",
                              "nevt total", "x"facs[1]" ["nevt1b[facs[1]]"]", "x"facs[2]" ["nevt1b[facs[2]]"]", "x"facs[3]" ["nevt1b[facs[3]]"]",
                              sabp1[tag], sap1[tag];
                              print lsep};
           printf "| %-10s |", tag;
           for(ifac=1; ifac<=nfac; ifac++)
           {fac=facs[ifac]; printf " %6.2f = %6.2f + %6.2f |", sec3[tag,fac,1], sec3[tag,fac,2], sec3[tag,fac,3]};
           printf " %8s |", tputm1[tag];
           if(tag=="FORTRAN"){ printf " %8s | %8s |", "---", "---"; }
           else{ printf " %8.2e |", tputb1[tag]; printf " %8.2e |", tput1[tag]; }
           printf "\n"};
          print lsep;
          print "";
         }' >> $out
  done
}

cd $(dirname $0)/..
echo PWD=$(pwd)

suff=mad
fpt=d
inl=inl0
hrd=hrd0

procs="eemumu ggtt ggttg ggttgg ggttggg"
###procs="ggttggg"

out=tmad/summaryTable_default.txt
oneTable
cat $out
