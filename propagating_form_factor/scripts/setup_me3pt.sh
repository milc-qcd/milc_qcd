#! /bin/sh

# Create input file for mat_el.x for extracting form-factor
# C. DeTar 6/16/98
# modified 12/31/98

# Usage (example)

#     setup_me3pt.sh <matelname> <p> <sq> <zk>

# Requires parameter files: param.<matelname> and q momentum list qlist


case $# in
[0-3]) echo "Usage $0 <matelname> <p> <sq> <zk>"
  exit 1
esac

matelname=$1
p3=$2
sq=$3
zk=$4

param=param.${matelname}

if [ ! -f ${param} ]
then
  echo "Missing parameter file ${param}"
  exit 1
fi

q3list=qlist

for q3 in `cat ${q3list}`
do
  inputparam=in.${matelname}_p${p3}_sq${sq}_zk${zk}_q${q3}

  . ${param}

  cat /dev/null > ${inputparam}

  echo "three_point_select ${fb3pt}" >> ${inputparam}
  echo "SP 0 ZK ${zk} SQ ${sq} Q ${q3} P 0 OP ${op3} CP 0 WT 1." >> ${inputparam}

# Look up correct recoil momentum in addition table

  k2=`awk '{if($1=="P" && $2 == '${p3}' && $3 == "Q" && $4 == '${q3}' && $5 == "K"){print $6;exit;}}' qpklist`

  echo "" >>  ${inputparam}
  echo "two_point_recoil_select ${fb2pt2}" >> ${inputparam}
  echo "SP 0 ZK ${zk} K ${k2} OP ${op2} CP 0 WT 1." >> ${inputparam}

  # Look up correct initial momentum index in table

  p1=`awk '{if($1=="P" && $2 == '${p3}' && $3 == "Q" && $4 == 0 && $5 == "K"){print $6;exit;}}' qpklist`

  echo "two_point_sequential_select ${fb2pt2}" >> ${inputparam}
  echo "SP 0 SQ ${sq} P ${p1} OP ${op1} CP 0 WT 1." >> ${inputparam}

  # Name of list file is derived from the 1st 6 letters of the file name
  list3=`echo $name3 | awk '{print substr($1,1,6)}'`
  list3=${list3}list

  echo "" >>  ${inputparam}
  echo "filelist ${list3}" >> ${inputparam}

done
