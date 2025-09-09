#!/bin/sh

# Run after running compile_quda.sh
# Run as:
#	bash compile_ks_spectrum_hisq.sh

################### Setup environment for Frontier

module reset
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
module load craype-accel-amd-gfx90a
module load cray-mpich/8.1.28
module load cmake
module load perftools
module load ninja
module list

QUDA_INSTALL=`pwd`/build/usqcd

export PATH=${ROCM_PATH}/bin:${ROCM_PATH}/llvm/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${ROCM_PATH}/llvm/lib64:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${ROCM_PATH}/llvm/lib:${MPICH_DIR}/lib:${GTL_ROOT}:${LD_LIBRARY_PATH}

# Can't seem to find mpi.h otherwise:
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/opt/cray/pe/mpich/8.1.28/ofi/amd/5.0/include

LIBQUDA="$(pat_opts pre_link hipcc gpu) -Wl,-rpath ${QUDA_INSTALL}/lib -L${QUDA_INSTALL}/lib -lquda -D__gfx90a --amdgpu-target=gfx90a -Wl,-rpath=${ROCM_PATH}/hiprand/lib -L${ROCM_PATH}/hiprand/lib -Wl,-rpath=${ROCM_PATH}/rocfft/lib -L${ROCM_PATH}/rocfft/lib -lhiprand -lrocfft -Wl,-rpath=${ROCM_PATH}/hipblas/lib -L${ROCM_PATH}/hipblas/lib -lhipblas -Wl,-rpath=${ROCM_PATH}/rocblas/lib -L${ROCM_PATH}/rocblas/lib -lrocblas -Wl,-rpath=${ROCM_PATH}/hip/lib $(pat_opts post_link hipcc gpu) -g -pg"

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
OPT="$(pat_opts include hipcc gpu) $(pat_opts pre_compile hipcc gpu)  -g -pg -ggdb -O3 -Ofast --offload-arch=gfx90a $(pat_opts post_compile hipcc gpu)" \
PATH_TO_NVHPCSDK="" \
CUDA_HOME="" \
LDFLAGS=" --verbose -L/opt/cray/pe/mpich/8.1.28/ofi/amd/5.0/lib -lmpi" \
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
