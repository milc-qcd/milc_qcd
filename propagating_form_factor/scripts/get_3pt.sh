#! /bin/sh

# Extracts raw three-pt functions
# CD 6/16/98
# modified 12/31/98

# Usage (example)

#     get_3pt.sh <matelname>  <p> <sp> <sq> <zk>

# Requires parameter files: param.<inputparam>

case $# in
[0-3])
  echo "Usage $0  <matelname>  <p> <sp> <sq> <zk>"
  exit 1
esac

matelname=$1
p3=$2
sp=$3
sq=$4
zk=$5

outpath=/work/u2219/prop_form/b560m01/3pt

param=param.${matelname}

. param.${matelname}

# Name of list file is derived from the 1st 6 letters of the file name
# Plus the spectator quark designation
list3stem=`echo $name3 | awk '{print substr($1,1,6)}'`
list3=${list3stem}_sp${sp}list

inputparam=in3.${matelname}_p${p3}_sp${sp}_sq${sq}_zk${zk}
threept_out=${outpath}/me.${matelname}_p${p3}_sp${sp}_sq${sq}_zk${zk}

for file in `cat ${list3}`
do
  # Process only files with names that match the stem
  if [ -n "`echo ${file} | grep ${list3stem}`" ]
  then
    cfg=`echo $file | awk -F. '{print $NF}'`
    multiselect_3pt.x  ${file} < ${inputparam}
  fi
done

gzip ${threept_out}*
/bin/rm ${inputparam}




