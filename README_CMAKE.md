# Use CMake to Configure MILC Project

## Build `su3_rhmd_hisq` with MPI and OpenMP

```bash
cd milc_qcd
mkdir -p build && cd build
cmake .. -DMPP=ON -DOMP=ON -DSUBDIR=ks_imp_rhmc -DPRECISION=1
cmake --build . -j8 --target su3_rhmd_hisq
ldd ks_imp_rhmc/su3_rhmd_hisq
```

## Build `su3_rhmd_hisq` with MPI, OpenMP and QUDA

```bash
cd milc_qcd
mkdir -p build && cd build
cmake .. -DMPP=ON -DOMP=ON -DSUBDIR=ks_imp_rhmc -DPRECISION=1 \
    -DWANTQUDA=ON \
    -DWANT_FN_CG_GPU=ON \
    -DWANT_FL_GPU=ON \
    -DWANT_GF_GPU=ON \
    -DWANT_FF_GPU=ON \
    -DWANT_GA_GPU=ON \
    -DQUDA_DIR=/path-to-quda-build/lib/cmake/lib/QUDA
cmake --build . -j8 --target su3_rhmd_hisq
ldd ks_imp_rhmc/su3_rhmd_hisq
```

## Build `ks_specturm_hisq` with MPI, OpenMP and QUDA

```bash
cd milc_qcd
mkdir -p build && cd build
cmake .. -DMPP=ON -DOMP=ON -DSUBDIR=ks_spectrum -DPRECISION=1 \
    -DWANTQUDA=ON \
    -DWANT_FN_CG_GPU=ON \
    -DWANT_FL_GPU=ON \
    -DWANT_GF_GPU=ON \
    -DWANT_FF_GPU=ON \
    -DWANT_GA_GPU=ON \
    -DQUDA_DIR=/path-to-quda-build/lib/cmake/lib/QUDA
cmake --build . -j8 --target su3_rhmd_hisq
ldd ks_specturm/ks_specturm_hisq
```
