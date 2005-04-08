#! /bin/sh

# Usage 

#  seterr.sh PAT1 PAT2 target tol


PAT1=$1
PAT2=$2
target=$3
tol=$4

outtest=out.test.$target
outsamp=out.sample.$target
outerrt=out.errtol.$target
outtesttmp=$outtest.tmp
outsamptmp=$outsamp.tmp

../headtail.pl $PAT1 $PAT2 < $outtest > $outtesttmp
../headtail.pl $PAT1 $PAT2 < $outsamp > $outsamptmp
../seterrfile.pl $outtesttmp $outsamptmp $tol > $outerrt

