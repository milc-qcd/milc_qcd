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
    -DGPU_FN_CG=ON \
    -DGPU_FL=ON \
    -DGPU_GF=ON \
    -DGPU_FF=ON \
    -DGPU_GA=ON \
    -DQUDA_DIR=/path-to-quda-build/lib/cmake/lib/QUDA
cmake --build . -j8 --target su3_rhmd_hisq
ldd ks_imp_rhmc/su3_rhmd_hisq
```

## Build `ks_spectrum_hisq` with MPI, OpenMP and QUDA

```bash
cd milc_qcd
mkdir -p build && cd build
cmake .. -DMPP=ON -DOMP=ON -DSUBDIR=ks_spectrum -DPRECISION=1 \
    -DWANTQUDA=ON \
    -DGPU_FN_CG=ON \
    -DGPU_FL=ON \
    -DGPU_GF=ON \
    -DGPU_FF=ON \
    -DGPU_GA=ON \
    -DQUDA_DIR=/path-to-quda-build/lib/cmake/lib/QUDA
cmake --build . -j8 --target ks_spectrum_hisq
ldd ks_spectrum/ks_spectrum_hisq
```

## Limitations

- QOP/QDP, Hadrons, and QPhiXJ dependency integration is not implemented.
  `WANTHADRONS` and `WANTQPHIXJ` are placeholders.
- Grid runtime settings are fixed in `CMakeLists.txt`, not cache options.
- Grid gauge-force sources are absent; use `GPU_GF=OFF` for CPU gauge force.
- Some legacy applications still reference missing sources or have compilation errors.
