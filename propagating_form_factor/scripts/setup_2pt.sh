#! /bin/sh

# Prepares input files for two-pt functions

# CD 6/16/98
# modified 1/8/98

# Usage (example)

#     setup_2pt.sh <twoptname> <sp> <zk>

# Requires parameter files: param.<twoptname> and k momentum list klist

# Label for matrix element
case $# in
[0-2]) echo "Usage $0 <twoptname> <zk>"
  exit 1
esac

twoptname=$1
sp=$2
zk=$3

# Momentum transfer list
klist=klist
# Where output ASCII files go
outpath=/work/u2219/prop_form/b560m01/2pt

inputparam=in2.${twoptname}_sp${sp}_zk${zk}
. param.${twoptname}

cat /dev/null > ${inputparam}

for k in `cat ${klist}`
do
  twopt_out=${outpath}/ms.${twoptname}_sp${sp}_zk${zk}_k${k}
  . in.${twoptname}
done
