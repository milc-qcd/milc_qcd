#! /bin/sh

# Prepares input files for raw three-pt functions - mainly for checking
# effective mass plateaus
# CD 6/16/98
# modified 1/8/98

# Usage (example)

#     setup_3pt.sh <matelname>  <p> <sp> <sq> <zk>

# Requires parameter files: param.<matelname> and q momentum list qlist


case $# in
[0-4]) echo "Usage $0 <matelname> <p> <sp> <sq> <zk>"
  exit 1
esac

matelname=$1
p3=$2
sp=$3
sq=$4
zk=$5

param=param.${matelname}
outpath=/work/u2219/prop_form/b560m01/3pt

if [ ! -f ${param} ]
then
  echo "Missing parameter file ${param}"
  exit 1
fi

q3list=qlist

inputparam=in3.${matelname}_p${p3}_sp${sp}_sq${sq}_zk${zk}
. ${param}

cat /dev/null > ${inputparam}

for q3 in `cat ${q3list}`
do
  threept_out=${outpath}/me.${matelname}_p${p3}_sp${sp}_sq${sq}_zk${zk}_q${q3}
  . in.${matelname}
done
