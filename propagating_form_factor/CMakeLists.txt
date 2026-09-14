# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    calc_heavy_light_form_factor3.c
    clover_rotate_fermilab.c
    clover_rotate.c
    contract_heavy_to_heavy.c
    contract_heavy_to_light.c
    contract_two_pt_prop.c
    copy_site_spin_wilson_vector.c
    corr_prop_routines.c
    debug_form.c
    dslash_lean_space.c
    gauge_info.c
    meson_spectrum.c
    sequential_source.c
    setup_control_form.c
    setup_form.c
    setup_timeslice_fft.c
    sink_smear_quark.c
    write_propagating_form_lean.c
    write_propagating_two_pt_lean.c
    zonked_quark_routines.c
    # TODO: Restore the source referenced by Make_template.
    d_w_meson.c
    # TODO: Restore the source referenced by Make_template.
    wilson_invert_lean.c
)

add_interface_library(g_objects
    ../generic_form/c_scale_wilson_vector2.c
    ../generic/layout_timeslices.c
    ../generic/bsd_sum.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic_wilson/io_helpers_w.c
    ../generic/io_lat4.c
    ../generic_wilson/io_prop_w.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/restrict_fourier.c
    ../generic_form/load_momentum.c
    ../generic_form/load_smearing.c
    ../generic_form/make_gammas.c
    ../generic_form/msw_gamma_l.c
    ../generic_form/msw_gamma_r.c
    ../generic_form/w_meson_finite_mom_lean3.c
    ../generic_form/wilson_source_routines.c
    ../generic_wilson/bj_to_weyl.c
    ../generic_wilson/copy_site_wilson_vector.c
    ../generic_wilson/dslash_w3.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
    ../generic_wilson/wilson_invert.c
)

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    machine_dep_io
    commpkg
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(w_mr_prop_form.x control_form.c
    DEFINES CGTIME
    EXTRA_SOURCES
    wilson_info.c
    # TODO: Restore the source referenced by Make_template.
    mrilu_w_or.c
)

add_milc_executable(cl_bicg_prop_form.x control_form.c
    DEFINES CLOVER CGTIME LARGE_MEMORY_SMEAR
    EXTRA_SOURCES
    clover_info.c
    ../generic_clover/f_mu_nu.c
    ../generic_clover/make_clov2.c
    ../generic_clover/d_bicgilu_cl.c
    ../generic_clover/d_hopilu_cl.c
)

