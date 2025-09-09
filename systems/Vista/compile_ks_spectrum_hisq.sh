#!/bin/sh

# Run after running compile_quda.sh
# Run as:
#	bash compile_ks_spectrum_hisq.sh

pushd .

module reset
module load gcc cuda
module list

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

export C_INCLUDE_PATH=$C_INCLUDE_PATH:/home1/apps/nvidia24/openmpi5/fftw3/3.3.10/include

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
LDFLAGS="-g -L/home1/apps/nvidia/Linux_aarch64/24.9/cuda/lib64 -lcudart -L/home1/apps/nvidia/Linux_aarch64/24.9/math_libs/lib64 -lcublas -lcufft -L/home1/apps/nvidia/Linux_aarch64/24.9/cuda/lib64/stubs -lcuda -lnvidia-ml" \
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
WANTPRIMME=false \
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

