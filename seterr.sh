#! /bin/sh

# Usage 

#  seterr.sh PAT1 PAT2 target tol

if [ $# -lt 4 ] 
then
  echo "Usage $0 PAT1 PAT2 target tol"
  exit
fi

PAT1=$1
PAT2=$2
target=$3
tol=$4

outtest=out.test.$target
outsamp=out.sample.$target
outerrt=out.errtol.$target
outtesttmp=$outtest.tmp
outsamptmp=$outsamp.tmp

perl ../headtail.pl $PAT1 $PAT2 < $outtest > $outtesttmp
perl ../headtail.pl $PAT1 $PAT2 < $outsamp > $outsamptmp
perl ../seterrfile.pl $outtesttmp $outsamptmp $tol > $outerrt

