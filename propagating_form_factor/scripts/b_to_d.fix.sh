#! /bin/sh

sp=2

for r in "" _rs _rz
do
  cd /work/u2219/prop_form/b560m01/3pt
  tar -cf me.bvd${r}_p0_sp${sp}.tar me.bv?d${r}_p0_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p0_sp${sp}.tar`
  tar -cf me.bvd${r}_p1_sp${sp}.tar me.bv?d${r}_p1_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p1_sp${sp}.tar`
  tar -cf me.bvd${r}_p2_sp${sp}.tar me.bv?d${r}_p2_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p2_sp${sp}.tar`
  cd $HOME/v5/propagating_form_factor/scripts
done
