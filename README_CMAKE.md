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

## Targets disabled because sources are missing

The following CMake target declarations are commented out, including targets
that inherit missing sources through shared groups. This inventory uses the
default CPU options; optional backend limitations above still apply. Paths are
relative to the application directory. The traditional Make templates retain their original
targets; restoring a source requires reviewing the corresponding CMake comments.

| Directory | Commented targets | Missing files |
| --- | --- | --- |
| `ks_imp_rhmc` | `test_su3_paranoid`, `test_su3_paranoid_der` | `control_test_su3_paranoid.c` |
| `ks_imp_utilities` | `check_ks_invert_hisq` | `control_ks_invert.c`, `setup_ks_invert.c` |

## Applications without CMake support

The following applications have no `CMakeLists.txt`. Upstream marked them
"unsuitable for use" (Carleton DeTar, 2011-11-29) or excluded them from the
supported list in `Make_test_all`. Except for `ks_dynamical`, they also do not
build with the traditional Make system because they reference removed sources.
`w_source.c` and `ks_source.c` were removed in 2011 in favor of the
`quark_source` API in `generic/quark_source*.c`; `layout_hyper.c` was removed
in 2005.

| Directory | Upstream status | Missing files |
| --- | --- | --- |
| `arb_dirac_eigen` | Commented out in `Make_test_all` | `../generic_wilson/w_source.c`, `layout_hyper.c` |
| `arb_dirac_invert` | Commented out in `Make_test_all` | `../generic_wilson/w_source.c`, `layout_hyper.c` |
| `clover_hybrids` | Commented out in `Make_test_all` | `../generic_wilson/w_source.c`, `ks_source.c` |
| `clover_invert` | Commented out in `Make_test_all` | `../generic_wilson/w_source.c`, `control_cl_hl_multi.c` |
| `dense_static_su3` | Not suitable for use (`bf46efc9`) | `layout_hyper.c`, `density_effective.c`, `external_field.c` |
| `h_dibaryon` | Not suitable for use (`e337ea6c`) | `../generic_wilson/w_source.c`, `dslash_lean.c`, `layout_hyper.c`, `make_clov.c`, `wilson_invert_lean.c` |
| `heavy` | Not suitable for use (`6e1587e1`) | `ks_source.c` |
| `hqet_heavy_to_light` | Unsuitable for use (`20be60d2`) | `../generic_wilson/w_source.c`, `dslash_lean.c`, `layout_hyper.c`, `make_clov.c`, `mrilu_w_or.c`, `wilson_invert_lean.c` |
| `ks_dynamical` | Unsuitable for use (`798d57f4`) | None |
| `ks_hl_spectrum` | Unsuitable for use (`8087531d`) | `../generic_wilson/w_source.c`, `ks_source.c` |
| `propagating_form_factor` | Unsuitable for use, not supported (`70779ccb`) | `../generic_wilson/w_source.c`, `d_w_meson.c`, `mrilu_w_or.c`, `wilson_invert_lean.c` |
| `pw_nr_meson` | Unsuitable for use (`d78be39a`) | `../generic_wilson/w_source.c`, `ks_source.c` |
| `schroed_pg` | Unsuitable for use (`a51332b2`) | `layout_hyper.c` |
| `string_break` | Unsuitable for use (`a51332b2`) | `layout_hyper.c` |
| `wilson_dynamical` | Unsuitable for use (`3dfc41a8`) | `d_congrad2.c`, `dslash_w.c`, `layout_hyper.c` |
