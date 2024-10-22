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

if [ "$MULTIGRID" = "1" ]
then
  MG="-DMULTIGRID"
fi

if [ ! -f "./Makefile" ]
then
  cp ../Makefile .
fi

if [ -f "./ks_spectrum_hisq" ]
then
  rm ./ks_spectrum_hisq
fi

#ARCH="pow9" \
#COMPILER="gnu" \
#OPT="-O3 -Ofast" \
#CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME" \
MY_CC=mpicc \
MY_CXX=mpicxx \
CUDA_HOME=${PATH_TO_CUDA} \
QUDA_HOME=${PATH_TO_QUDA} \
WANTQUDA=true \
WANT_FN_CG_GPU=true \
WANT_FL_GPU=true \
WANT_GF_GPU=true \
WANT_FF_GPU=true \
WANT_GAUGEFIX_OVR_GPU=true \
WANT_MIXED_PRECISION_GPU=2 \
PRECISION=2 \
MPP=true \
OMP=true \
WANTQIO=true \
WANTQMP=true \
QIOPAR=${PATH_TO_QIO} \
QMPPAR=${PATH_TO_QMP} \
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID -DMULTISOURCE $MG" \
make -j 1 ks_spectrum_hisq

