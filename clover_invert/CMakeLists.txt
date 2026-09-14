# Source groups and targets follow Make_template.

add_interface_library(my_objects
    clover_info.c
    gauge_info.c
    setup_cl.c
    spectrum_cl.c
    ../generic_wilson/dirac_utilities.c
)

add_interface_library(g_objects
    ../generic/ape_smear.c
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/discretize_wf.c
    ../generic/field_utilities.c
    ../generic/field_translation.c
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
    ../generic/ranstuff.c
    ../generic/remap_stdio_from_args.c
    ../generic/report_invert_status.c
    ../generic/reunitarize2.c
    ../generic_wilson/baryon_cont.c
    ../generic_wilson/canopy2weyl_rot.c
    ../generic_wilson/dslash_w_space.c
    ../generic_wilson/dslash_w3.c
    ../generic_wilson/gammas.c
    ../generic_wilson/gauss_smear_w.c
    ../generic_wilson/io_helpers_w.c
    ../generic_wilson/io_prop_w.c
    ../generic_wilson/io_prop_w_fm.c
    ../generic_wilson/io_source_w_fm.c
    ../generic_wilson/meson_cont.c
    ../generic_wilson/staggered2naive.c
    ../generic_wilson/w_baryon.c
    ../generic_wilson/wilson_invert.c
    ../generic_wilson/w_meson.c
    # TODO: Restore the source referenced by Make_template.
    ../generic_wilson/w_source.c
    ../generic_clover/f_mu_nu.c
    ../generic_clover/make_clov2.c
    ../generic/layout_hyper_prime.c
    ../generic_clover/cl_solver_utilities.c
    ../generic_clover/d_cgilu_cl.c
    ../generic_clover/d_hopilu_cl.c
    ../generic_clover/d_mrilu_cl.c
    ../generic_clover/d_bicgilu_cl.c
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
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
        ../generic_ks/io_scidac_ks.c
        ../generic_wilson/io_scidac_w.c
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

add_milc_executable(su3_clov_bi control_cl.c
    DEFINES IOTIME PRTIME
)
target_compile_options(su3_clov_bi PRIVATE -UCL_CG -DCL_CG=BICG)

add_milc_executable(su3_clov_hop control_cl.c
    DEFINES IOTIME PRTIME
)
target_compile_options(su3_clov_hop PRIVATE -UCL_CG -DCL_CG=HOP)

add_milc_executable(su3_clov_mr control_cl.c
    DEFINES IOTIME PRTIME
)
target_compile_options(su3_clov_mr PRIVATE -UCL_CG -DCL_CG=MR)

add_milc_executable(su3_clov_cg control_cl.c
    DEFINES IOTIME PRTIME
)
target_compile_options(su3_clov_cg PRIVATE -UCL_CG -DCL_CG=CG)

add_milc_executable(su3_clhl_cg control_cl_hl.c
    DEFINES HL
    EXTRA_SOURCES
    spectrum_cl_hl.c
    ../generic_wilson/w_baryon_hl.c
)
target_compile_options(su3_clhl_cg PRIVATE -UCL_CG -DCL_CG=CG)

add_milc_executable(su3_clhl_bi control_cl_hl.c
    DEFINES HL
    EXTRA_SOURCES
    spectrum_cl_hl.c
    ../generic_wilson/w_baryon_hl.c
)
target_compile_options(su3_clhl_bi PRIVATE -UCL_CG -DCL_CG=BICG)

add_milc_executable(su3_clhl_hop control_cl_hl.c
    DEFINES HL
    EXTRA_SOURCES
    spectrum_cl_hl.c
    ../generic_wilson/w_baryon_hl.c
)
target_compile_options(su3_clhl_hop PRIVATE -UCL_CG -DCL_CG=HOP)

add_milc_executable(su3_clhl_hop_multi control_cl_hl_multi.c
    DEFINES HL
    EXTRA_SOURCES
    spectrum_cl_hl.c
    ../generic_wilson/w_baryon_hl.c
)
target_compile_options(su3_clhl_hop_multi PRIVATE -UCL_CG -DCL_CG=HOP)

