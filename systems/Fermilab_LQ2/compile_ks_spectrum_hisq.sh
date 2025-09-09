#!/bin/sh

# Run after running compile_quda.sh
# Run as:
#	bash compile_ks_spectrum_hisq.sh

pushd .

module reset
module load gompi cuda ucc_cuda ucx_cuda
module list

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Trick to automatically determine CUDA_HOME
CUDA_PATH_TMP=`which nvcc`
CUDA_HOME=${CUDA_PATH_TMP/\/bin\/nvcc/}

# QUDA and USQCD build directories: Modify paths as needed
cd build
QUDA_BUILD=`pwd`
USQCD_BUILD=${QUDA_BUILD}/usqcd
cd ..

# Download/update MILC
if [ -d milc_qcd ]
then
  cd milc_qcd/milc_qcd
  git checkout develop
  git pull
else
  mkdir milc_qcd
  cd milc_qcd
  git clone https://github.com/milc-qcd/milc_qcd.git
  cd milc_qcd
  git checkout develop
fi

############ Make ks_spectrum_hisq ##################
cd ks_spectrum
cp ../Makefile .
make clean

MY_CC=mpicc \
MY_CXX=mpiCC \
ARCH="" \
GPU_ARCH="nvidia" \
OFFLOAD="CUDA" \
COMPILER="gnu" \
OPT="-O3 -Ofast -g" \
LDFLAGS="-g -L/srv/software/el8/x86_64/hpc/cuda/12.2.1/lib64/stubs/ -lcuda -lnvidia-ml " \
CUDA_HOME=${CUDA_HOME} \
QUDA_HOME=${QUDA_BUILD} \
WANTQUDA=true \
WANT_FN_CG_GPU=true \
WANT_FL_GPU=true \
WANT_GF_GPU=true \
WANT_FF_GPU=true \
WANT_MIXED_PRECISION_GPU=2 \
PRECISION=2 \
WANT_GAUGEFIX_OVR_GPU=true \
WANT_GSMEAR_GPU=true \
MPP=true \
OMP=true \
WANTQIO=true \
WANTQMP=true \
QIOPAR=${USQCD_BUILD} \
QMPPAR=${USQCD_BUILD} \
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID " \
CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DFLTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME -DGS_TIME" \
make -j 1 ks_spectrum_hisq >& make_ks_spectrum_hisq.log

popd

echo ""
echo ""
echo "Check that compilation was successful by viewing make_ks_baryon.log:"
echo ""
tail milc_qcd/milc_qcd/ks_spectrum/make_ks_spectrum_hisq.log

