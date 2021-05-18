#!/bin/bash

if [ `uname -p` != 'x86_64' ] ; then
    echo "host architecture is not x86_64!"
    return
fi

module load gnu8 openmpi3 cuda10/10.2

export CUDA_HOME=/srv/software/cuda-toolkits/10.2
export CUBLAS_LIBDIR=/srv/software/cuda-toolkits/nvidia-hpc/Linux_x86_64/20.7/math_libs/10.2/lib64
echo QUDA_HOME=${QUDA_HOME}

MPP=false
OMP=false

if [ ${MPP} == true ] ; then
    module load openmpi3
    export QUDA_HOME=/work1/simone/contract/quda_c7340d4_x86_64
    MY_CC=mpicc
    MY_CXX=mpicxx
    WANTQMP=true
    QMPPAR=${QUDA_HOME}/usqcd
    QIOPAR=${QUDA_HOME}/usqcd
    QMPSNG=unused
    QIOSNG=unused
else
    export QUDA_HOME=/work1/simone/contract/quda_c7340d4_serial_x86_64
    MY_CC=gcc
    MY_CXX=g++
    WANTQMP=false
    QMPSNG=/work1/simone/contract/scidac_x86_64/install/qmp-single
    QIOSNG=/work1/simone/contract/scidac_x86_64/install/qio-single
    QMPPAR=unused
    QIOPAR=unused
fi
module list

MY_CC=${MY_CC} \
    MY_CXX=${MY_CXX} \
    CUDA_HOME=${CUDA_HOME} \
    QUDA_HOME=${QUDA_HOME} \
    CUBLAS_LIBDIR=${CUBLAS_LIBDIR} \
    WANTQUDA=true \
    WANT_FN_CG_GPU=true \
    WANT_FL_GPU=true \
    WANT_GF_GPU=true \
    WANT_FF_GPU=true \
    WANT_MIXED_PRECISION_GPU=2 \
    WANT_KS_CONT_GPU=true \
    WANT_KS_CONT_CPU=false \
    PRECISION=2 \
    MPP=${MPP} \
    OMP=${OMP} \
    WANTQIO=true \
    WANTQMP=${WANTQMP} \
    QIOPAR=${QIOPAR} \
    QIOSNG=${QIOSNG} \
    QMPPAR=${QMPPAR} \
    QMPSNG=${QMPSNG} \
    make -f Makefile.x86_64 -j16 ks_spectrum_hisq

