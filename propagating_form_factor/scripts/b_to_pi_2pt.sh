#! /bin/sh

sp=2
for zk in 0 1 2 3 4 
do
  setup_2pt.sh b_GG ${sp} ${zk}
  get_2pt.sh   b_GG ${sp} ${zk}
  setup_2pt.sh b_GE ${sp} ${zk}
  get_2pt.sh   b_GE ${sp} ${zk}
  setup_2pt.sh b_GL ${sp} ${zk}
  get_2pt.sh   b_GL ${sp} ${zk}
  setup_2pt.sh p_GG ${sp} ${zk}
  get_2pt.sh   p_GG ${sp} ${zk}
  setup_2pt.sh r_GG ${sp} ${zk}
  get_2pt.sh   r_GG ${sp} ${zk}
  cd /work/u2219/prop_form/b560m01/2pt
  tar -cf ms.zk${zk}.tar  ms.b_GG_sp${sp}_zk${zk}_k??.gz ms.b_GE_sp${sp}_zk${zk}_k??.gz ms.b_GL_sp${sp}_zk${zk}_k??.gz ms.p_GG_sp${sp}_zk${zk}_k??.gz ms.r_GG_sp${sp}_zk${zk}_k??.gz
  rm `tar -tf ms.zk${zk}.tar`
  cd $HOME/v5/propagating_form_factor/scripts
done
