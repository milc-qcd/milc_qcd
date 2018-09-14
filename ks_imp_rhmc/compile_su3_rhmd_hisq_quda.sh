#!/bin/bash

CC=mpicc \
CXX=mpicxx \
CUDA_HOME=/usr/local/cuda-9.2 \
QUDA_HOME=/scratch/kate/quda-develop \
WANTQUDA=true \
WANT_FN_CG_GPU=true \
WANT_FL_GPU=true \
WANT_GF_GPU=true \
WANT_FF_GPU=true \
PRECISION=2 \
MPP=true \
OMP=true \
make -j 1 su3_rhmd_hisq
