#! /bin/sh

sp=0

for r in "" _rs _rz
do
  for p in 0 1 2
  do
    for sq in 0 1 2 3 4
    do
      for d in x y z t
      do
	setup_3pt.sh bv${d}d${r} ${p} ${sp} ${sq} ${sq}
	get_3pt.sh   bv${d}d${r} ${p} ${sp} ${sq} ${sq}
      done
    done
  done

  cd /work/u2219/prop_form/b560m01/3pt
  tar -cf me.bvd${r}_p0_sp${sp}.tar me.bv?d${r}_p0_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p0_sp${sp}.tar`
  tar -cf me.bvd${r}_p1_sp${sp}.tar me.bv?d${r}_p1_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p1_sp${sp}.tar`
  tar -cf me.bvd${r}_p2_sp${sp}.tar me.bv?d${r}_p2_sp${sp}_sq*
  rm `tar -tf me.bvd${r}_p2_sp${sp}.tar`
  cd $HOME/v5/propagating_form_factor/scripts
done
