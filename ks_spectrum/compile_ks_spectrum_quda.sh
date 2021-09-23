
if [ -z "$PATH_TO_CUDA" ]
then
  echo "Environment variable PATH_TO_CUDA unset, exiting..."
fi

if [ -z "$PATH_TO_QUDA" ]
then
  echo "Environment variable PATH_TO_QUDA unset, exiting..."
fi

if [ -z "$PATH_TO_QIO" ]
then
  echo "Environment variable PATH_TO_QIO unset, exiting..."
fi

if [ -z "$PATH_TO_QMP" ]
then
  echo "Environment variable PATH_TO_QMP unset, exiting..."
fi

if [ ! -f "./Makefile" ]
then
  cp ../Makefile .
fi

if [ -f "./ks_spectrum_hisq" ]
then
  rm ./ks_spectrum_hisq
fi

# Remove CTIME to remove overly verbose timing output. Useful for debugging.
# Update verbosity
sed -i 's/CGPU += -DSET_QUDA_SUMMARIZE/CGPU += -DSET_QUDA_VERBOSE/g' Makefile

if [[ ! -z "$POWER9" ]] && [[ "$POWER9" -eq "1" ]]
then
  export ARCH="pow9"
  export COMPILER="gnu"
  export OPT="-O3 -Ofast"
fi

#ARCH="pow9" \
#COMPILER="gnu" \
#OPT="-O3 -Ofast" \
MY_CC=mpicc \
MY_CXX=mpicxx \
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
CGEOM="-DFIX_NODE_GEOM -DFIX_IONODE_GEOM" \
KSCGMULTI="-DKS_MULTICG=HYBRID -DMULTISOURCE -DMULTIGRID" \
CTIME="-DNERSC_TIME -DCGTIME -DFFTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME" \
make -j 1 ks_spectrum_hisq

