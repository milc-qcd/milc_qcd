#! /bin/sh

SCRATCH=/work/u2219/prop_form/b560m01

for matelname in bvtp
do
  r=rs; p=1 ; sq=1
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rs; p=1 ; sq=2
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz; p=0 ; sq=4
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz;
  for p in 1 2
  do
    for sq in 0 1 2 3 4
    do
      rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
      setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
      get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0
    done
  done

done

for matelname in bvxp
do
  r=rs; p=1 ; sq=3
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz;
  for p in 1 2
  do
    for sq in 0 1 2 3 4
    do
      rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
      setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
      get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0
    done
  done

done

for matelname in bvyp
do
  r=rs; p=1 ; sq=3
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz; p=0 ; sq=4
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz;
  for p in 1 2
  do
    for sq in 0 1 2 3 4
    do
      rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
      setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
      get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0
    done
  done

done


for matelname in bvzp
do
  r=rz; p=0 ; sq=4
  rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
  get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0

  r=rz;
  for p in 1 2
  do
    for sq in 0 1 2 3 4
    do
      rm $SCRATCH/me.${matelname}_${r}_p${p}_sq${sq}_zk0_q*
      setup_3pt.sh ${matelname}_${r} ${p} ${sq} 0
      get_3pt.sh   ${matelname}_${r} ${p} ${sq} 0
    done
  done

done

for matelname in bvtp bvyp bvzp
do
  p=1 ; sq=0
  rm $SCRATCH/me.${matelname}_p${p}_sq${sq}_zk0_q*
  setup_3pt.sh ${matelname} ${p} ${sq} 0
  get_3pt.sh   ${matelname} ${p} ${sq} 0

done

for matelname in bvxp
do
  p=1
  for sq in 1 2
  do
    rm $SCRATCH/me.${matelname}_p${p}_sq${sq}_zk0_q*
    setup_3pt.sh ${matelname} ${p} ${sq} 0
    get_3pt.sh   ${matelname} ${p} ${sq} 0
  done

done

