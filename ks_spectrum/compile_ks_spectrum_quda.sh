
PATH_TO_CUDA=/sw/summit/cuda/10.1.243
QUDA_BUILD=/gpfs/alpine/proj-shared/lgt104/detar/kpp_spectrum/l192384f211b728m000415m01129m1329-mg/new_runs/build-optimized
USQCD_BUILD=${QUDA_BUILD}/usqcd

# Uncomment ARCH, COMPILER, OPT for power9
# Remove CTIME to remove overly verbose timing output. Useful for debugging.

# Update verbosity
sed -i 's/CGPU += -DSET_QUDA_SUMMARIZE/CGPU += -DSET_QUDA_VERBOSE/g' Makefile

ARCH="pow9" \
COMPILER="gnu" \
OPT="-O3 -Ofast" \
MY_CC=mpicc \
MY_CXX=mpicxx \
CUDA_HOME=${PATH_TO_CUDA} \
QUDA_HOME=${QUDA_BUILD} \
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
QIOPAR=${USQCD_BUILD} \
QMPPAR=${USQCD_BUILD} \
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID -DMULTISOURCE -DMULTIGRID" \
CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME" \
make -j 1 ks_spectrum_hisq

