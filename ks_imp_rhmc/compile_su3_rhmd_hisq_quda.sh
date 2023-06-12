#!/bin/bash

if [ -z "$PATH_TO_CUDA" ]
then
  echo "Environment variable PATH_TO_CUDA unset, exiting..."
  exit
fi

if [ -z "$PATH_TO_QUDA" ]
then
  echo "Environment variable PATH_TO_QUDA unset, exiting..."
  exit
fi

if [ -z "$PATH_TO_QIO" ]
then
  echo "Environment variable PATH_TO_QIO unset, exiting..."
  exit
fi

if [ -z "$PATH_TO_QMP" ]
then
  echo "Environment variable PATH_TO_QMP unset, exiting..."
  exit
fi

if [ ! -f "./Makefile" ]
then
  cp ../Makefile .
fi

if [ -f "./su3_rhmd_hisq" ]
then
  rm ./su3_rhmd_hisq
fi

#ARCH="pow9" \
#COMPILER="gnu" \
#OPT="-O3 -Ofast" \
#CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DGATIME -DGFTIME -DREMAP -DPRTIME -DIOTIME" \
#QUDA_VERBOSITY=VERBOSE \
MY_CC=mpicc \
MY_CXX=mpicxx \
CUDA_HOME=${PATH_TO_CUDA} \
QUDA_HOME=${PATH_TO_QUDA} \
WANTQUDA=true \
WANT_FN_CG_GPU=true \
WANT_FL_GPU=true \
WANT_GF_GPU=true \
WANT_FF_GPU=true \
WANT_GA_GPU=true \
WANT_MIXED_PRECISION_GPU=2 \
PRECISION=2 \
MPP=true \
OMP=true \
WANTQIO=true \
WANTQMP=true \
QIOPAR=${PATH_TO_QIO} \
QMPPAR=${PATH_TO_QMP} \
make -j 1 su3_rhmd_hisq
