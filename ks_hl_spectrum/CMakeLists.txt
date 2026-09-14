# Source groups and targets follow Make_template.

add_interface_library(my_objects
    clover_info.c
    gauge_info.c
    ksprop_info.c
    r_smear.c
    rotate_wq.c
    # TODO: Restore the source referenced by Make_template.
    ks_source.c
)

add_interface_library(g_objects
    ../generic_wilson/canopy2weyl_rot.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/discretize_wf.c
    ../generic/field_translation.c
    ../generic/field_utilities.c
    ../generic/gaugefix2.c
    ../generic/gauge_utilities.c
    ../generic/io_detect.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/io_source_cmplx_fm.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/ranstuff.c
    ../generic/remap_stdio_from_args.c
    ../generic/reunitarize2.c
    ../generic/restrict_fourier.c
    ../generic/ploop3.c
    ../generic/ape_smear.c
    ../generic_wilson/dirac_utilities.c
    ../generic_wilson/dslash_w_space.c
    ../generic_wilson/gammas.c
    ../generic_wilson/gauss_smear_w.c
    ../generic_wilson/io_helpers_w.c
    ../generic_wilson/io_prop_w.c
    ../generic_wilson/io_prop_w_fm.c
    ../generic_wilson/io_source_w_fm.c
    ../generic/project_su3_hit.c
    ../generic_wilson/staggered2naive.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
    ../generic_wilson/wilson_invert.c
    ../generic_ks/d_congrad5.c
    ../generic_ks/d_congrad_opt.c
    ../generic_ks/dslash.c
    ../generic_ks/flavor_ops2.c
    ../generic_ks/gauss_smear_ks.c
    ../generic_ks/io_helpers_ks.c
    ../generic_ks/io_prop_ks.c
    ../generic_ks/io_prop_ks_fm.c
    ../generic_ks/ks_utilities.c
    ../generic_ks/rephase.c
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

add_milc_executable(ks_hl_spectrum control_hl_spectrum.c
    EXTRA_SOURCES
    twopt.c
    setup.c
)

add_milc_executable(ks_hl_baryon control_hl_baryon.c
    EXTRA_SOURCES
    setup_baryon.c
    baryon_twopt.c
)

