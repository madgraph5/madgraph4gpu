#!/bin/sh
# Copyright (C) 2020-2024 CERN and UCLouvain.
# Licensed under the GNU Lesser General Public License (version 3 or later).
# Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
# Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

scrdir=$(cd $(dirname $0); pwd)

echo "==============================================================================="
echo "Executing $0 $*"
echo "==============================================================================="
echo

# Include CUDA/8tpb?
###cuda8tpb=
cuda8tpb="CUDA/8tpb"

# Quiet?
quiet=
if [ "$1" == "-quiet" ]; then
  quiet=$1; shift
fi

# Short tables?
onlyxmax=1
if [ "$1" == "-long" ] && [ "$2" != "-ALL" ]; then
  onlyxmax=0
  shift
fi

# Choose table(s)
table=
if [ "$1" == "-ALL" ] && [ "$2" == "" ]; then
  set -e
  ###$0 ${quiet} -long -default
  ###$0 ${quiet} -default
  $0 ${quiet} -v10000
  exit 0
###elif [ "$1" == "-default" ]; then
###  table="default"; shift
elif [ "$1" == "-v10000" ]; then
  table="v10000"; shift
else
  echo "Usage: $0 [-quiet] [-long] <table [-ALL|-default|-v10000]>"; exit 1
fi

# Select revisions and characteristics of mad logs
# Select processes and fptypes
mrevs=""
###if [ "$table" == "default" ]; then
###  procs=...
###  mrevs=...
###  fpts=...
###  taglist=...
if [ "$table" == "v10000" ]; then
  #procs="ggttgg ggttggg"
  procs="ggttgg"
  taglist="FORTRAN CPP/none CPP/sse4 CPP/avx2 CPP/512y CPP/512z CUDA/8192 CUDA/max $cuda8tpb"
  mrevs="$mrevs cd8e872" # cuda120/gcc113  (03 Oct 2024 itscrd90) v1.00.00
  #mrevs="$mrevs a3d64bd" # -------/gcc114  (03 Oct 2024 gold91)   v1.00.00
  #mrevs="$mrevs 07c2a53" # roc60/clang170  (03 Oct 2024 lumi)     v1.00.00+amd
  fpts="d f m"
fi
revs="$mrevs"

# TEST MODE (debug individual lines and skip the final END printout)
testmode=0
###procs="ggttggg"; testmode=1

# Kernel function
function oneTable()
{
  parlist="(1) (2-none) (2-sse4) (2-avx2) (2-512y) (2-512z) (3-cuda) (3bis)"
  faclist="1 10"
  for proc in $procs; do
    file=tmad/logs_${proc}_${suff}/log_${proc}_${suff}_${fpt}_${inl}_${hrd}.txt
    if [ ! -f $file ]; then continue; fi
    ###echo "*** FILE $file ***"
    git checkout $rev $file >& /dev/null
    if [ "$?" != "0" ]; then echo "ERROR! 'git checkout $rev' failed!"; exit 1; fi
    node=$(cat $file | grep ^On | sort -u)
    if [ "$nodelast" != "$node" ]; then echo -e "$node\n" >> $out; nodelast=$node; fi
    cat $file | awk -vproc=$proc -vtaglist="$taglist" -vparlist="$parlist" -vfaclist="$faclist" -vonlyxmax=$onlyxmax -vtestmode=$testmode '
      BEGIN{status=testmode} # set status=1 to skip the END printout and debug individual lines 
      BEGIN{ntag=split(taglist,tags); npar=split(parlist,pars);
            ###if(ntag!=npar){print "ERROR! ntag!=npar", ntag, npar; status=1; exit status}; # NB new ntag>npar!
            for(i=1;i<=npar;i++){tag1[pars[i]]=tags[i];}}
      BEGIN{nfac=split(faclist,facs)}
      BEGIN{if(onlyxmax==0) lsepEQUAL=sprintf("%0136d",0); else lsepEQUAL=sprintf("%0107d",0);
            lsepDASH=lsepEQUAL; gsub("0","-",lsepDASH); gsub("0","=",lsepEQUAL)}
      BEGIN{if(onlyxmax==0) lsepEQUAL2=sprintf("%014d%97s%025d",0,"",0); else lsepEQUAL2=sprintf("%014d%68s%025d",0,"",0);
            lsepDASH2=lsepEQUAL2; gsub("0","-",lsepDASH2); gsub("0","=",lsepEQUAL2)}
      BEGIN{if(onlyxmax==0) ifac0=1; else ifac0=nfac}
      ###/create events.lhe/{print $0}
      /create events.lhe/{par=$2; tag=tag1[par]} # current tag (FORTRAN... CUDA/8192)
      ###/create events.lhe/{print par, tag}
      /GCHECK\(MAX\)/{tag="CUDA/max"} # current tag (CUDA/max)
      /GCHECK\(MAX128THR\)/{tag="CUDA/max128t"} # current tag (CUDA/max128t)
      /GCHECK\(MAX8THR\)/{tag="CUDA/8tpb"} # current tag (CUDA/8tpb)
      /create events.lhe/{fac=substr($5,2)} # current fac
      /\[XSECTION\] nb_page_loop =/{if(tag!="") nloop2[tag,fac]=$4}
      ###/\[COUNTERS\]/{print $0}
      /\[COUNTERS\]/{typadd=0; # NEW (TEMPORARY?) add HEL to Overhead
                     if($3=="TOTAL") typ=1;
                     else if($3=="Overhead") typ=2;
                     else if($3=="MEs") typ=3;
                     else if($3=="HEL") {typ=4; typadd=2} # NEW (TEMPORARY?) add HEL to Overhead
                     else{print "ERROR! Unknown type "$3; status=1; exit status};
                     if($4==":") sec=$5; else sec=$8; sec=substr(sec,1,length(sec)-1); # current sec
                     sec3[tag,fac,typ]=sec;
                     if(typadd>0) sec3[tag,fac,typadd]+=sec} # NEW (TEMPORARY?) add HEL to Overhead
      /\[COUNTERS\]/{if(tag!="" && $3=="MEs")
                     {nevt=$10; ###print tag, nevt;
                      if(tag=="FORTRAN") nevt1[fac]=nevt; 
                      else if(nevt1[fac]!=nevt){print "ERROR! nevt mismatch", nevt1[fac], nevt; status=1; exit status};
                      tputm1[tag]=tolower($15); # TODO? cross-check this against nevt1[fac]/sec3[tag,fac,3]
                      tputm1tot[tag]=nevt/sec3[tag,fac,1]}}
      /CHECK/{if(tag!=""){if($5=="2048") fld5="2k";
                          else if($5=="16384") fld5="16k";
                          else if($5=="65536") fld5="64k";
                          else fld5=$5;
                          if($8=="--bridge"){gcheck="bridge"; sabg1[tag]=$5*$6; sabp1[tag]=fld5"*"$6"*"$7}
                          else{gcheck="nobridge"; sag1[tag]=$5*$6; sap1[tag]=fld5"*"$6"*"$7}}}
      ###/CHECK/{print $0}
      ###/CHECK/{print tag, sabg1[tag], sabp1[tag], sag1[tag], sap[tag]; }
      /EvtsPerSec/{if(tag!="" && gcheck=="bridge"){tputb1[tag]=$5}}
      /EvtsPerSec/{if(tag!="" && gcheck!="bridge"){tput1[tag]=$5}}
      END{if (status!=0) exit status;
          print lsepEQUAL;
          if(onlyxmax==0)
            printf "| %-10s | mad%23s | mad%23s | mad%14s | mad%14s | %-9s | %-9s |\n",
                   "", "x"facs[1], "x"facs[2], "x"facs[2], "x"facs[2], "sa/brdg", "sa/full";
          else
            printf "| %-10s | mad%23s | mad%14s | mad%14s | %-9s | %-9s |\n",
                   ###"", "x"facs[2], "x"facs[2], "x"facs[2], "sa/brdg", "sa/full";
                   "", "", "", "", "sa/brdg", "sa/full";
          print lsepDASH;
          if(onlyxmax==0)
            printf "| %-10s | %-26s | %-26s | %-17s | %-17s | %-9s | %-9s |\n",
                   proc, "[sec] tot = mad + MEs", "[sec] tot = mad + MEs",
                   "[TOT/sec]", "[MEs/sec]", "[MEs/sec]", "[MEs/sec]";
          else
            printf "| %-10s | %-26s | %-17s | %-17s | %-9s | %-9s |\n",
                   proc, "[sec] tot = mad + MEs",
                   "[TOT/sec]", "[MEs/sec]", "[MEs/sec]", "[MEs/sec]";
          print lsepEQUAL;
          for (itag=1; itag<=ntag; itag++)
          {tag=tags[itag]; 
           if(tag=="FORTRAN"){if(onlyxmax==0)
                                printf "| %-10s | %26s | %26s | %17s | %17s | %9s | %9s |\n",
                                "nevt/grid", nloop2[tag,fac], nloop2[tag,fac], nloop2[tag,fac], nloop2[tag,fac],
                                sabg1["CPP/none"], sag1["CPP/none"];
                              else
                                printf "| %-10s | %26s | %17s | %17s | %9s | %9s |\n",
                                "nevt/grid", nloop2[tag,fac], nloop2[tag,fac], nloop2[tag,fac],
                                sabg1["CPP/none"], sag1["CPP/none"];
                              if(onlyxmax==0)
                                printf "| %-10s | %26s | %26s | %17s | %17s | %9s | %9s |\n",
                                "nevt total", nevt1[facs[1]], nevt1[facs[2]], nevt1[facs[2]], nevt1[facs[2]],
                                sabp1["CPP/none"], sap1["CPP/none"];
                              else
                                printf "| %-10s | %26s | %17s | %17s | %9s | %9s |\n",
                                "nevt total", nevt1[facs[2]], nevt1[facs[2]], nevt1[facs[2]],
                                sabp1["CPP/none"], sap1["CPP/none"];
                              print lsepDASH}
           else if(tag=="CUDA/max"||tag=="CUDA/8tpb"){
                              if(tag=="CUDA/max") print lsepEQUAL; else print lsepEQUAL2;
                              if(onlyxmax==0)
                                printf "| %-10s | %95s | %9s | %9s |\n",
                                "nevt/grid", "", sabg1[tag], sag1[tag];
                              else
                                printf "| %-10s | %66s | %9s | %9s |\n",
                                "nevt/grid", "", sabg1[tag], sag1[tag];
                              if(onlyxmax==0)
                                printf "| %-10s | %95s | %9s | %9s |\n",
                                "nevt total", "", sabp1[tag], sap1[tag];
                              else
                                printf "| %-10s | %66s | %9s | %9s |\n",
                                "nevt total", "", sabp1[tag], sap1[tag];
                              print lsepDASH2};
           printf "| %-10s |", tag;
           if(tag=="CUDA/max"||tag=="CUDA/8tpb")
                            { if(onlyxmax==0) printf " %95s |", "";
                              else printf " %66s |", ""; }
           else{ for(ifac=ifac0; ifac<=nfac; ifac++)
                 { fac=facs[ifac]; printf " %7.2f = %6.2f + %7.2f |", sec3[tag,fac,1], sec3[tag,fac,2], sec3[tag,fac,3]};
                 if(tag=="FORTRAN") txttot="="; else txttot="x";
                 ###if(tag=="CPP/none") txtmes="="; else txtmes="x";
                 if(tag=="FORTRAN") txtmes="="; else txtmes="x";
                 ratiotot=sprintf("%4.1f",tputm1tot[tag]/tputm1tot["FORTRAN"]);
                 if(length(ratiotot)>4) ratiotot=substr(ratiotot,0,4);
                 ###ratiomes=sprintf("%4.1f",tputm1[tag]/tputm1["CPP/none"]);
                 ratiomes=sprintf("%4.1f",tputm1[tag]/tputm1["FORTRAN"]);
                 if(length(ratiomes)>4) ratiomes=substr(ratiomes,0,4);
                 printf " %9.2e (%1s%4s) | %9s (%1s%4s) |", 
                        tputm1tot[tag], txttot, ratiotot,
                        tputm1[tag], txtmes, ratiomes; }
           if(tag=="FORTRAN"){ printf " %9s | %9s |", "---", "---"; }
           else{ printf " %9.2e |", tputb1[tag]; printf " %9.2e |", tput1[tag]; }
           if(tag=="CUDA/max"||tag=="CUDA/8tpb")
                            { printf "\n| %-10s |", "";
                              if(onlyxmax==0) printf " %95s |", "";
                              else printf " %66s |", "";
                              ratiomes2=sprintf("%4.1f",tput1[tag]/tputm1["FORTRAN"]);
                              if(length(ratiomes2)>4) ratiomes2=substr(ratiomes2,0,4);
                              printf " %-9s |", ""; printf "   (x%4s) |", ratiomes2; }
           printf "\n"};
	  if(tag=="CUDA/max"||tag=="CUDA/8tpb") print lsepEQUAL2; else print lsepEQUAL;
          print "\n";
         }' >> $out
  done
}

cd $(dirname $0)/..
###echo PWD=$(pwd)

suff=mad
inl=inl0
hrd=hrd0

if [ "${onlyxmax}" == "1" ]; then
  out=${scrdir}/summaryTable_${table}.txt
else
  out=${scrdir}/summaryTable_${table}_long.txt
fi  
echo "" > $out
for fpt in $fpts; do
  echo -e "*** FPTYPE=$fpt ******************************************************************\n" >> $out
  for rev in $revs; do
    echo -e "+++ $bckend REVISION $rev (commit date: $(git log $rev --pretty=format:'%ci' --abbrev-commit -n1)) +++" >> $out
    oneTable
  done
done

git checkout HEAD ../cudacpp/tmad/logs* >& /dev/null
if [ "$quiet" == "" ]; then cat $out; fi
