#! /bin/sh

sp=2
r="_rs"
for zk in 2
do
  for p in 0 1 2
  do
    for sq in 0 1 2 3 4
    do
      setup_3pt.sh bvtp${r} ${p} ${sp} ${sq} ${zk}
      get_3pt.sh   bvtp${r} ${p} ${sp} ${sq} ${zk}
      setup_3pt.sh bvxp${r} ${p} ${sp} ${sq} ${zk}
      get_3pt.sh   bvxp${r} ${p} ${sp} ${sq} ${zk}
      setup_3pt.sh bvyp${r} ${p} ${sp} ${sq} ${zk}
      get_3pt.sh   bvyp${r} ${p} ${sp} ${sq} ${zk}
      setup_3pt.sh bvzp${r} ${p} ${sp} ${sq} ${zk}
      get_3pt.sh   bvzp${r} ${p} ${sp} ${sq} ${zk}
    done
  done
  cd /work/u2219/prop_form/b560m01/3pt
  tar -cf me.bvp${r}_p0_sp${sp}_zk${zk}.tar me.bv?p${r}_p0_sp${sp}_sq?_zk${zk}*
  rm `tar -tf me.bvp${r}_p0_sp${sp}_zk${zk}.tar`
  tar -cf me.bvp${r}_p1_sp${sp}_zk${zk}.tar me.bv?p${r}_p1_sp${sp}_sq?_zk${zk}*
  rm `tar -tf me.bvp${r}_p1_sp${sp}_zk${zk}.tar`
  tar -cf me.bvp${r}_p2_sp${sp}_zk${zk}.tar me.bv?p${r}_p2_sp${sp}_sq?_zk${zk}*
  rm `tar -tf me.bvp${r}_p2_sp${sp}_zk${zk}.tar`
  cd $HOME/v5/propagating_form_factor/scripts
done

