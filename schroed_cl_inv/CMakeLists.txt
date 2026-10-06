# Source groups and targets follow Make_template.

add_interface_library(my_objects
    clover_info.c
    setup_cl.c
    schroed_meson.c
    ape_smear_SF.c
    zv_meas.c
    gauge_info.c
)

add_interface_library(g_objects
    ../generic/check_unitarity.c
    ../generic/d_plaq4.c
    ../generic/field_utilities.c
    ../generic/io_helpers.c
    ../generic/io_lat4.c
    ../generic/io_lat_utils.c
    ../generic/make_lattice.c
    ../generic/nersc_cksum.c
    ../generic/plaquette4.c
    ../generic/remap_stdio_from_args.c
    ../generic/report_invert_status.c
    ../generic/reunitarize2.c
    ../generic_wilson/dslash_lean2.c
    ../generic_wilson/gauss_smear_w.c
    ../generic_wilson/io_source_w_fm.c
    ../generic_wilson/meson_cont.c
    ../generic_wilson/wilson_invert.c
    ../generic_wilson/canopy2weyl_rot.c
    ../generic_clover/f_mu_nu.c
    ../generic_clover/make_clov2.c
    ../generic_schroed/make_schroed_lattice.c
    ../generic/layout_hyper_prime.c
    ../generic_clover/cl_solver_utilities.c
    ../generic_clover/d_cgilu_cl.c
    ../generic_clover/d_hopilu_cl.c
    ../generic_clover/d_mrilu_cl.c
    ../generic_clover/d_bicgilu_cl.c
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
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

add_milc_executable(su3_schr_cl_cg control_cl.c
    DEFINES
)
target_compile_options(su3_schr_cl_cg PRIVATE -UCL_CG -DCL_CG=CG)

add_milc_executable(su3_schr_cl_bi control_cl.c
    DEFINES
)
target_compile_options(su3_schr_cl_bi PRIVATE -UCL_CG -DCL_CG=BICG)

add_milc_executable(su3_schr_cl_mr control_cl.c
    DEFINES
)
target_compile_options(su3_schr_cl_mr PRIVATE -UCL_CG -DCL_CG=MR)

