# Source groups and targets follow Make_template.
# Legacy machine-dependent I/O uses the current serial/MPI communication
# and ANSI I/O interfaces; old source names are retained below.

add_interface_library(my_objects
    3submat.c
    apply_hqet_proj.c
    calc_hqet_light_form_factor.c
    clover_rotate.c
    contract_hqet_to_light.c
    contract_two_pt.c
    copy_spin_wilson_vector.c
    corr_routines.c
    debug_form.c
    evolve_hqet_prop.c
    gauge_info.c
    hqet_prop.c
    hqet_prop_backwards.c
    load_velocity.c
    m_mat_an_z.c
    m_mat_nn_z_inc.c
    sequential_hqet_source.c
    setup_control_hqet_form.c
    setup_hqet_form.c
    setup_timeslice_fft.c
    smear_hqet_prop.c
    twopt_sequential.c
    write_hl_hqet_form.c
    write_seq_twopt.c
    write_twopt.c
    # TODO: Restore the source referenced by Make_template.
    layout_hyper.c
    # TODO: Restore the source referenced by Make_template.
    dslash_lean.c
    # TODO: Restore the source referenced by Make_template.
    wilson_invert_lean.c
)

add_interface_library(g_objects
    ../generic_form/c_scale_wilson_vector2.c
    ../generic/bsd_sum.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic_wilson/io_helpers_w.c
    ../generic/io_lat4.c
    ../generic_wilson/io_prop_w.c
    ../generic/make_lattice.c
    ../generic/restrict_fourier.c
    ../generic_form/flip_source_rep.c
    ../generic_form/load_momentum.c
    ../generic_form/load_smearing.c
    ../generic_form/make_gammas.c
    ../generic_form/msw_gamma_l.c
    ../generic_form/msw_gamma_r.c
    ../generic_form/w_meson_finite_mom.c
    ../generic_wilson/bj_to_weyl.c
    ../generic_wilson/copy_site_wilson_vector.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
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

add_milc_executable(w_mr_hqet_form.x control_form.c
    EXTRA_SOURCES
    wilson_info.c
    # TODO: Restore the source referenced by Make_template.
    mrilu_w_or.c
)

add_milc_executable(cl_bicg_hqet_form.x control_form.c
    DEFINES BICG_CLOVER
    EXTRA_SOURCES
    clover_info.c
    ../generic_clover/f_mu_nu.c
    # TODO: Restore the source referenced by Make_template.
    make_clov.c
    ../generic_clover/d_bicgilu_cl.c
)

