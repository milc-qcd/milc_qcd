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
- There is no LAPACK option (Make uses a hand-set `LIBLAPACK`), so the eigCG
  targets cannot link unless LAPACK comes in through another package.
- `ARCH=native` is a CMake-only choice and adds no compiler flags.
- Some legacy applications still reference missing sources or have compilation
  errors; see [Known issues inherited from Make](#known-issues-inherited-from-make).

## Known issues inherited from Make

The CMake build mirrors `Makefile`, `Make_template_combos`, and each
application's `Make_template`, including their mistakes. The issues below
are known and intentionally left unchanged so that both build systems produce
the same programs. The build results below use the default CPU options
(`cmake -DSUBDIR=<dir>`).

### Build selection

- `WANTQUDA=ON` with `GPU_FN_CG=OFF` fails to compile. With QUDA, `KSCGSTORE`
  drops `DBLSTORE_FN`, but the CPU branch still selects
  `generic_ks/dslash_fn_dblstore.c`, which requires it.
- The multimass inverter selects QPhiX on `HAVE_QPHIX` (package enabled), while
  the single-mass CG selects it on `HAVE_FN_CG_QPHIX` (`WANT_FN_CG_QPHIX`). With
  `WANTQPHIX=ON` and `WANT_FN_CG_QPHIX=OFF`, single-mass CG uses MILC and
  multimass CG uses QPhiX.
- With Grid and `GPU_FF=ON`, only `hisq_force` is defined; asqtad targets have
  no `asq_force`.
- `ARCH=epyc` only affects `GRID_ARCH` in Make, so it adds no compiler flags.
- `crc32` and `get_qmp_node_number` are always defined, but they need the QIO
  and QMP headers (`WANTQIO`, `WANTQMP`).

### Targets that need optional packages

| Directory | Targets | Requires |
| --- | --- | --- |
| `ext_src` | `ext_src` | `WANTQIO` |
| `file_combine` | `su3_combine` | `WANTQIO` |
| `file_utilities` | `crc32`, `diff_wprop`, `lattice_to_scidac`, `make_modulation`, `scidac_2v5`, `v5_2scidac` | `WANTQIO` |
| `file_utilities` | `get_qmp_node_number` | `WANTQMP` |
| `ks_imp_utilities` | `check_reunit_deriv` | `WANTQIO` |
| `ks_imp_utilities` | `make_links_hisq` | `WANTQIO` or `WANTGRID` |
| `ks_measure` | `ks_measure_current_hisq`, `ks_measure_current_hisq_u1`, `ks_measure_current_hisq_u1_loop`, `ks_measure_current_eigcg_hisq` | `WANTQIO` |
| `ks_measure` | `ks_measure_eigcg_hisq` | LAPACK |
| `ks_spectrum` | `ks_spectrum_eigcg_asqtad`, `ks_spectrum_eigcg_hisq` | LAPACK |
| `u1_gauge` | `u1_g`, `u1_g_convert` | `WANTFFTW` |

### Targets that fail because of source or Make template errors

| Directory | Targets | Cause |
| --- | --- | --- |
| `arb_overlap` | `su3_ov_eig_cg_f_hyp`, `su3_ov_eig_cg_f_hyp_per`, `su3_ov_eig_cg_multi`, `su3_ov_eig_cg_multi_per` | `arb_ov_includes.h:176`: `#else HAVE_ARPACK` followed by a second `#else` |
| `clover_dynamical` | `su3_hmc_bi_screen`, `su3_hmc_spectrum`, `su3_rmd_spectrum` | Missing comma after `F_OFFSET(psi)` in `wilson_invert_site` calls (`s_props_cl.c:61`, `w_spectrum_cl.c:64`) |
| `clover_dynamical` | `su3_phi` | Defines `BI` but lists no `d_bicgilu_cl` object; `bicgilu_cl_site` is undefined |
| `file_utilities` | `ckpt1_to_v5` | Lists no generic I/O objects; `big_endian`, `complete_U`, and others are undefined |
| `file_utilities` | `diff_colormatrix`, `diff_ksprop`, `diff_ksvector` | Includes KS headers but sets no `QUARK`, so `quark_action.h` is missing |
| `gluon_prop` | `su3_asqtad_quark_prop`, `su3_fn_quark_prop`, `su3_ks_quark_prop` | `generic_ks/mat_invert.c` uses `eigVec`, `eigVal`, and `param`, which the application does not declare |
| `gluon_prop` | `su3_asqtad_renorm` | `quark_renorm.c:220`: `dslash_site(...` is not closed |
| `gluon_prop` | `su3_gluon_prop_imp` | `gaugefixfft.c:170`: `ftmp1` is undeclared |
| `gluon_prop` | `su3_p4_quark_prop` | `generic_ks/d_congrad5_eo.c`: prototypes conflict with the current headers |
| `hvy_qpot` | `su3_hqp_ape`, `su3_hqp_coulomb`, `su3_hqp_coulomb_plcor`, `su3_hqp_hyp`, `su3_hybrid_hqp_ape`, `su3_hybrid_hqp_hyp` | `generic/io_helpers.c:760` uses `site_prn`, which `lattice.h` does not declare |
| `ks_eigen` | `su3_eigen_p4fat3` | `generic_ks/fermion_links_milc.c`: `eo_links_t` has no `preserve` member |
| `ks_imp_dyn` | `su3_hmc_eo_symzk1_p4`, `su3_rmd_eo_symzk1_p4` | `generic_ks/fermion_links_milc.c`: `eo_links_t` has no `preserve` member |
| `ks_imp_dyn` | `su3_rmd_eo_symzk1_fat7tad` | `generic_ks/d_congrad5_eo.c`: prototypes conflict with the current headers |
| `ks_imp_rhmc` | `su3_rhmc_hisq_debug`, `su3_rhmc_hisq_su3_debug`, `su3_rhmc_hisq_tune`, `su3_rhmc_hisq_wrap_asq` | `update_h_rhmc.c:101` uses `fn_links.hl`; `fermion_links_t` is a pointer with no `hl` member |
| `ks_imp_rhmc` | `test_su3_mat_op_inverse` | `hisq/hisq_action.h` does not define `UNITARIZATION_GROUP` |
| `ks_imp_rhmc` | `test_su3_mat_op_unit`, `test_su3_mat_op_unit_der` | Calls `su3_unitarize_analytic` and `su3_unit_der_analytic`, which are declared only in some configurations |
| `ks_imp_utilities` | `check_fermion_force_asqtad`, `check_fermion_force_hisq` | `check_fermion_force.c:12` is a deliberate `BOMB` line |
| `ks_imp_utilities` | `check_link_fattening` | Make defines `CHECK_FATTENING`, but `setup.c` tests `CHECK_LINK_FATTENING` |
| `pure_gauge`, `pure_gauge_old` | `su3_ora_glue` | `generic/glueball_op.c:333`: `path_product` is called with 2 of 3 arguments |
| `rcorr` | `current_rcorr` | `lattice.h` uses `uint32_t` and `size_t` without including `<stdint.h>` |
| `schroed_cl_inv` | `su3_schr_cl_bi`, `su3_schr_cl_cg`, `su3_schr_cl_mr` | `generic/io_helpers.c:760` uses `site_prn`, which `lattice.h` does not declare |
| `schroed_ks_dyn` | `su3_schr_hmc`, `su3_schr_hmc_ph`, `su3_schr_rw_hmc` | `update.c:152`: missing `;` after `destroy_fn_links(fn)` |
| `schroed_ks_dyn` | `su3_schr_phi`, `su3_schr_phi_ph`, `su3_schr_rmd`, `su3_schr_rmd_ph`, `su3_schr_rw_phi`, `su3_schr_rw_rmd` | `generic_ks/mat_invert.c` uses `eigVec`, `eigVal`, and `param`, which the application does not declare |
| `schroed_ks_dyn` | `su3_schr_sw_hmc`, `su3_schr_sw_phi` | `update_hgh.c:92`: `get_fm_links` is called with 1 of 2 arguments |
| `smooth_inst` | `su3_ape`, `su3_fat`, `su3_hyp`, `su3_hyp2`, `su3_stout` | `save_topo.c` uses `int32type`, renamed to `u_int32type` upstream |

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
