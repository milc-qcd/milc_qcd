# Source groups and targets follow Make_template.

add_interface_library(my_objects
    overlap_info.c
    gauge_info.c
    readinfo.c
    setup_p_cl.c
    delta0y.c
    eigen_stuff_helpers.c
    eigen_stuff_PRIMME.c
    eigen_stuff_Ritz.c
    make_clov_field.c
    mult_ldu2.c
    f_mu_nu1.c
    path.c
    setup_offset.c
    setup_links.c
    hdelta0_field.c
    build_params_0.166.c
    wp_shrink_pl_field_pf.c
    wp_shrink_pl_field_pf_l.c
    wp_grow_pl_field_pf.c
    wp_grow_pl_field_pf_l.c
    grsource_dirac.c
    io_modes.c
    copy_fields.c
    myclock.c
    block_nhyp.c
    nhyp.c
    congrad_multi_field.c
    hoverlap.c
    build_hr0_field2.c
    build_h0.c
    build_lowest_chi.c
    build_hov.c
    vectorh.c
    jacobi.c
)

add_interface_library(g_objects
    ../generic/ploop3.c
    ../generic/ape_smear.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
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
    ../generic/momentum_twist.c
    ../generic/nersc_cksum.c
    ../generic/quark_source.c
    ../generic/quark_source_sink_op.c
    ../generic/quark_source_io.c
    ../generic/ranstuff.c
    ../generic/project_su3_hit.c
    ../generic/remap_stdio_from_args.c
    ../generic/report_invert_status.c
    ../generic/reunitarize2.c
    ../generic_wilson/boundary_flip.c
    ../generic_wilson/canopy2weyl_rot.c
    ../generic_wilson/dirac_utilities.c
    ../generic/discretize_wf.c
    ../generic_wilson/gammas.c
    ../generic_wilson/io_helpers_w.c
    ../generic_wilson/io_prop_w.c
    ../generic_wilson/io_prop_w_fm.c
    ../generic_wilson/io_source_w_fm.c
    ../generic_wilson/staggered2naive.c
    ../generic/layout_hyper_prime.c
)

if(HAVEFFTW)
    append_interface_library(g_objects
        ../generic/remap_fftw_fourier.c
    )
else()
    append_interface_library(g_objects
        ../generic/restrict_fourier.c
    )
endif()

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
        ../generic_wilson/io_scidac_w.c
        ../generic/file_types_milc_usqcd.c
    )
endif()

add_interface_library(objects
    OBJECTS
    my_objects
    g_objects
    cl_objects
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

target_compile_definitions(objects INTERFACE NHYP)

if(HAVE_PRIMME)
    target_compile_definitions(objects INTERFACE PRIMME)
endif()

add_milc_executable(su3_ov_eig_cg_f_hyp control_f.c
    DEFINES NHYP HAVE_DIRAC EIG EIGO FIELD RANDOM PROJ
)

add_milc_executable(su3_ov_eig_cg_multi control_f.c
    DEFINES NHYP HAVE_DIRAC EIG EIGO FIELD RANDOM PROJ INV MINN MULTI
    EXTRA_SOURCES
    multi_cg.c
    multi_cg_iter.c
)

add_milc_executable(su3_ov_eig_cg_f_hyp_per control_f.c
    DEFINES NHYP HAVE_DIRAC EIG EIGO FIELD RANDOM PROJ TPERIODIC
)

add_milc_executable(su3_ov_eig_cg_multi_per control_f.c
    DEFINES NHYP HAVE_DIRAC EIG EIGO FIELD RANDOM PROJ INV MINN MULTI TPERIODIC
    EXTRA_SOURCES
    multi_cg.c
    multi_cg_iter.c
)

