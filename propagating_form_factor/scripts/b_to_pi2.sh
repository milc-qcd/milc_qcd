#! /bin/sh

path=$HOME/v5/propagating_form_factor/scripts

#sp=1
#for zk in 0 1 2 3 4 
#do
#  setup_2pt.sh b_GG ${sp} ${zk}
#  get_2pt.sh   b_GG ${sp} ${zk}
#  setup_2pt.sh b_GE ${sp} ${zk}
#  get_2pt.sh   b_GE ${sp} ${zk}
#  setup_2pt.sh b_GL ${sp} ${zk}
#  get_2pt.sh   b_GL ${sp} ${zk}
#  setup_2pt.sh p_GG ${sp} ${zk}
#  get_2pt.sh   p_GG ${sp} ${zk}
#  setup_2pt.sh r_GG ${sp} ${zk}
#  get_2pt.sh   r_GG ${sp} ${zk}
#done

sp=1
#for r in "" _rs _rz
for r in _rz
do
#  for zk in 0 1 2
  for zk in 2
  do
#    for p in 0 1 2
    for p in 2
    do
#      for sq in 0 1 2 3 4
      for sq in 4
      do
#	 setup_3pt.sh bvtp${r} ${p} ${sp} ${sq} ${zk}
#	 get_3pt.sh   bvtp${r} ${p} ${sp} ${sq} ${zk}
#	 setup_3pt.sh bvxp${r} ${p} ${sp} ${sq} ${zk}
#	 get_3pt.sh   bvxp${r} ${p} ${sp} ${sq} ${zk}
#	 setup_3pt.sh bvyp${r} ${p} ${sp} ${sq} ${zk}
#	 get_3pt.sh   bvyp${r} ${p} ${sp} ${sq} ${zk}
	 setup_3pt.sh bvzp${r} ${p} ${sp} ${sq} ${zk}
	 get_3pt.sh   bvzp${r} ${p} ${sp} ${sq} ${zk}
      done
    done

    cd /work/u2219/prop_form/b560m01
    tar -cf met${r}${sp}${zk}.tar me.bvtp${r}_p*sp${sp}*zk${zk}*
    tar -cf mex${r}${sp}${zk}.tar me.bvxp${r}_p*sp${sp}*zk${zk}*
    tar -cf mey${r}${sp}${zk}.tar me.bvyp${r}_p*sp${sp}*zk${zk}*
    tar -cf mez${r}${sp}${zk}.tar me.bvzp${r}_p*sp${sp}*zk${zk}*
    cd $HOME/v5/propagating_form_factor/scripts

  done
done
