#! /bin/sh

# Compute the form-factor ratios
# C. DeTar 6/16/98
# modified 1/13/98

# Usage

#     get_me3pt.sh <inputparam>

if [ $# -lt 1 ]
then
  echo "Usage $0 <inputparam>"
  exit 1
fi

inputparam=$1

q3list=qlist

for q3 in `cat ${q3list}`
do
   ffout=`echo $inputparam | awk -F. '{print "ff." $2}'`

   # "input_param" is hard-wired in mat_el.x at the moment

   cp ${inputparam} input_param

   mat_el.x
   mv matrix_out ${ffout}
done
