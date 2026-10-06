#!/bin/sh

# Run after running compile_quda.sh
# Run as:
#	bash compile_ks_spectrum_hisq.sh

################### Setup environment for Frontier

module reset
module load PrgEnv-amd amd/7.1.1 rocm/7.1.1
module load craype-accel-amd-gfx90a
module load cmake
module load ninja
module list

QUDA_INSTALL=`pwd`/build/usqcd

MY_LDFLAGS="--verbose -g -Wl,-rpath=${MPICH_DIR}/lib -L${MPICH_DIR}/lib -lmpi"
MY_CFLAGS="-I${MPICH_DIR}/include -g -pg -ggdb -O3 -Ofast --offload-arch=gfx90a"

LIBQUDA="-Wl,-rpath ${QUDA_INSTALL}/lib -L${QUDA_INSTALL}/lib -lquda -D__gfx90a --amdgpu-target=gfx90a -Wl,-rpath=${ROCM_PATH}/hiprand/lib -L${ROCM_PATH}/hiprand/lib -Wl,-rpath=${ROCM_PATH}/rocfft/lib -L${ROCM_PATH}/rocfft/lib -lhiprand -lrocfft -Wl,-rpath=${ROCM_PATH}/hipblas/lib -L${ROCM_PATH}/hipblas/lib -lhipblas -Wl,-rpath=${ROCM_PATH}/rocblas/lib -L${ROCM_PATH}/rocblas/lib -lrocblas -Wl,-rpath=${ROCM_PATH}/hip/lib -g -pg"

echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "LIBQUDA: $LIBQUDA"

################### Get MILC

pushd .

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

################### Compile

cd ks_spectrum
cp ../Makefile .
make clean

OFFLOAD=HIP \
MY_CC=hipcc \
MY_CXX=hipcc \
COMPILER="gnu" \
ARCH="" \
OPT="${MY_CFLAGS}" \
PATH_TO_NVHPCSDK="" \
CUDA_HOME="" \
LDFLAGS=${MY_LDFLAGS} \
QUDA_HOME=${QUDA_INSTALL} \
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
QIOPAR=${QUDA_INSTALL} \
QMPPAR=${QUDA_INSTALL} \
LIBQUDA=${LIBQUDA} \
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID " \
CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DFLTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME -DGS_TIME" \
make -j 1 ks_spectrum_hisq >& make_ks_spectrum_hisq.log

popd

echo ""
echo ""
echo "Check that compilation was successful by viewing make_*.log:"
echo ""
tail milc_qcd/milc_qcd/ks_spectrum/make_ks_spectrum_hisq.log
