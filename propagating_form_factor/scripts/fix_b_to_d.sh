#! /bin/sh

for p in 1 2
do
for sq in 0 1 2 3 4
do
  setup_3pt.sh bvxd ${p} ${sq} 0
  get_3pt.sh   bvxd ${p} ${sq} 0
  setup_3pt.sh bvyd ${p} ${sq} 0
  get_3pt.sh   bvyd ${p} ${sq} 0
  setup_3pt.sh bvzd ${p} ${sq} 0
  get_3pt.sh   bvzd ${p} ${sq} 0
done
done


