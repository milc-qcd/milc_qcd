#! /bin/sh

# Extracts two-pt functions needed to reduce a three-pt function

# CD 6/16/98
# modified 12/31/98

# Usage (example)

#     get_2pt.sh <twoptname> <sp> <zk>

# Requires parameter files: param.<twoptname>

# Label for matrix element
case $# in
[0-2])
  echo "Usage $0 <twoptname> <zk>"
  exit 1
esac

twoptname=$1
sp=$2
zk=$3

outpath=/work/u2219/prop_form/b560m01/2pt

. param.${twoptname}

twoptlist=${name}_sp${sp}list

inputparam=in2.${twoptname}_sp${sp}_zk${zk}
meson_out=ms.${twoptname}_sp${sp}_zk${zk}

for file in `cat ${twoptlist}`
do
  multiselect_2pt.x  ${file} < ${inputparam}
done

gzip ${outpath}/${meson_out}*
/bin/rm ${inputparam}




