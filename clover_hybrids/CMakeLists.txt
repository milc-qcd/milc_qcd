# Source groups and targets follow Make_template.

add_interface_library(my_objects
    copy_site_wilson_vector.c
    f_measure_cl.c
    gauge_info.c
    ksprop_info.c
    mat_invert_clover.c
    setup.c
    smear_links_space.c
    spectrum_hybrids5_cl.c
    spectrum_light_light.c
    # TODO: Restore the source referenced by Make_template.
    ks_source.c
)

add_interface_library(g_objects
    ../generic/ape_smear.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_strength.c
    ../generic/gaugefix2.c
    ../generic/gauge_utilities.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/io_source_cmplx_fm.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/project_su3_hit.c
    ../generic/remap_stdio_from_args.c
    ../generic/ranstuff.c
    ../generic/reunitarize2.c
    ../generic_ks/io_helpers_ks.c
    ../generic_ks/io_prop_ks.c
    ../generic_ks/io_prop_ks_fm.c
    ../generic_ks/ks_utilities.c
    ../generic_wilson/baryon_cont.c
    ../generic_wilson/boundary_flip.c
    ../generic_wilson/canopy2weyl_rot.c
    ../generic_wilson/dirac_utilities.c
    ../generic_wilson/gauss_smear_w.c
    ../generic_wilson/io_source_w_fm.c
    ../generic_wilson/meson_cont.c
    ../generic_wilson/w_baryon.c
    ../generic_wilson/w_meson.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
    ../generic_wilson/wilson_invert.c
    ../generic_clover/f_mu_nu.c
    ../generic_clover/make_clov2.c
    ../generic/layout_hyper_prime.c
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
        ../generic_wilson/io_scidac_w.c
        ../generic_ks/io_scidac_ks.c
    )
endif()

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

add_milc_executable(hyb_bicg_clover control.c
    DEFINES LU BI SMEAR=1
    EXTRA_SOURCES
    ../generic_wilson/dslash_w3.c
    ../generic_clover/d_bicgilu_cl.c
    ../generic_clover/cl_solver_utilities.c
)

add_milc_executable(hyb_cg_clover_for control.c
    DEFINES LU SMEAR=1 FORWARDS_ONLY
    EXTRA_SOURCES
    ../generic_wilson/dslash_lean_forwards.c
    ../generic_clover/d_cgilu_cl.c
    ../generic_clover/cl_solver_utilities.c
)

