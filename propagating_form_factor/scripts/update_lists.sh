#! /bin/sh
#  Updates lists from the master list in HL3_p0_sp?list and HL2_GG_sp?list

for sp in 2
do
  sed 's/p0/p1/' HL3_p0_sp${sp}list > HL3_p1_sp${sp}list
  sed 's/p0/p2/' HL3_p0_sp${sp}list > HL3_p2_sp${sp}list
  for p in 0 1 2
  do
    sed 's/HL3/HH3/' HL3_p${p}_sp${sp}list | sed 's/LL2/HL2/' > HH3_p${p}_sp${sp}list
  done

  sed 's/GG/GE/' HL2_GG_sp${sp}list > HL2_GE_sp${sp}list
  sed 's/GG/GL/' HL2_GG_sp${sp}list > HL2_GL_sp${sp}list
  sed 's/HL/LL/' HL2_GG_sp${sp}list > LL2_GG_sp${sp}list
done
