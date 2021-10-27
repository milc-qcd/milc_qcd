#!/bin/bash

CC=mpicc \
CXX=mpicxx \
CUDA_HOME=${PATH_TO_CUDA} \
QUDA_HOME=${PATH_TO_QUDA} \
WANTQUDA=true \
WANT_FN_CG_GPU=true \
WANT_FL_GPU=true \
WANT_GF_GPU=true \
WANT_FF_GPU=true \
WANT_MIXED_PRECISION_GPU=2 \
PRECISION=2 \
MPP=true \
OMP=true \
WANTQIO=true \
WANTQMP=true \
QIOPAR=${PATH_TO_QIO} \
QMPPAR=${PATH_TO_QMP} \
make -j 1 su3_rhmd_hisq
