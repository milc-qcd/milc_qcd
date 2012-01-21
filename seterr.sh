#! /bin/sh

# Usage 

#  seterr.sh target tol PAT1 PAT2 PAT1 PAT2 ...

if [ $# -lt 4 ] 
then
  echo "Usage $0 target tol PAT1 PAT2 PAT1 PAT2 ..."
  exit
fi

target=$1
shift
tol=$1
shift

outtest=$target.test-out
outsamp=$target.sample-out
outerrt=$target.errtol
outtesttmp=$target.test-tmp
outsamptmp=$target.sample-tmp

perl ../../headtail.pl $* < $outtest > $outtesttmp
perl ../../headtail.pl $* < $outsamp > $outsamptmp
perl ../../seterrfile.pl $outtesttmp $outsamptmp $tol > $outerrt

