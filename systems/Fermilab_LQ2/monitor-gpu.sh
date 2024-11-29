#!/bin/bash

while [ true ]
do
  echo "monitor-gpu.sh reporting from `hostname`"
  date
  nvidia-smi
  sleep 300
done
