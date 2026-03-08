add_interface_library(my_objects
)

add_interface_library(g_objects
    ../generic/io_detect.c
    ../generic/io_lat_utils.c
    ../generic/remap_stdio_from_args.c
    ../generic/field_utilities.c
)

add_interface_library(g_ks_objects
)

add_interface_library(g_w_objects
)

if(HAVEQIO)
    append_interface_library(g_objects
        ../generic/file_types_milc_usqcd.c
        ../generic/io_scidac.c
        ../generic/io_scidac_types.c
    )
    append_interface_library(g_ks_objects
        ../generic_ks/io_scidac_ks.c
    )
    append_interface_library(g_w_objects
        ../generic_wilson/io_scidac_w.c
    )
endif()

add_interface_library(layout
    ../generic/layout_hyper_prime.c
)

add_interface_library(add_objects
    OBJECTS
    machine_dep_io
    commpkg
)

add_interface_library(objects
    OBJECTS
    my_objects
    layout
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(check_gauge check_gauge.c
    DEFINES ONLY_GLUON_FILES
    EXTRA_OBJECTS g_objects
    EXTRA_SOURCES
    ../generic/nersc_cksum.c
    ../generic/gauge_utilities.c
    ../generic/check_unitarity.c
    gauge_info_dummy.c
    ../generic/d_plaq4.c
    ../generic/io_lat4.c
)

add_milc_executable(diff_colormatrix diff_colormatrix.c
    EXTRA_OBJECTS g_objects g_ks_objects
    EXTRA_SOURCES
    gauge_info_dummy.c
    ../generic/gauge_utilities.c
    ../generic/d_plaq4.c
    ../generic/io_lat4.c
    ../generic/nersc_cksum.c

)

add_milc_executable(diff_ksprop diff_ksprop.c
    DEFINES NO_GAUGE_FIELD HAVE_KS
    EXTRA_OBJECTS g_objects g_ks_objects
    EXTRA_SOURCES
    ../generic/discretize_wf.c
    ../generic_wilson/gammas.c
    gauge_info_dummy.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic_ks/io_helpers_ks.c
    ../generic_ks/io_prop_ks.c
    ../generic_ks/io_prop_ks_fm.c
    ../generic/io_source_cmplx_fm.c
    ../generic_ks/io_source_ks_fm.c
    ksprop_info.c
    ../generic/quark_source.c
    ../generic/quark_source_io.c
    ../generic/quark_source_sink_op.c
    ../generic/ranstuff.c
    ../generic/restrict_fourier.c
    ../generic_ks/shift_field.c
    ../generic_ks/spin_taste_ops.c
)

add_milc_executable(diff_ksvector diff_ksvector.c
    EXTRA_OBJECTS g_objects g_ks_objects
    EXTRA_SOURCES
    ksprop_info.c
    ../generic/gauge_utilities.c
    gauge_info_dummy.c
    ../generic/d_plaq4.c
    ../generic/io_lat4.c
    ../generic/nersc_cksum.c
)

add_milc_executable(diff_wprop diff_wprop.c
    DEFINES NO_GAUGE_FIELD HAVE_DIRAC
    EXTRA_OBJECTS g_objects g_ks_objects g_w_objects
    EXTRA_SOURCES
    ../generic_wilson/canopy2weyl_rot.c
    ../generic/discretize_wf.c
    ../generic_wilson/gammas.c
    ../generic_ks/flavor_ops2.c
    gauge_info_dummy.c
    ../generic_wilson/gauss_smear_w.c
    ../generic/gauge_utilities.c
    ../generic/io_helpers.c
    ../generic_ks/io_helpers_ks.c
    ../generic_wilson/io_helpers_w.c
    ../generic_ks/io_prop_ks_fm.c
    ../generic_ks/io_prop_ks.c
    ../generic_wilson/io_prop_w.c
    ../generic_wilson/io_prop_w_fm.c
    ../generic/io_source_cmplx_fm.c
    ../generic_ks/io_source_ks_fm.c
    ../generic_wilson/io_source_w_fm.c
    ksprop_info.c
    ../generic/project_su3_hit.c
    ../generic/quark_source.c
    ../generic/quark_source_io.c
    ../generic/quark_source_sink_op.c
    ../generic/ranstuff.c
    ../generic/restrict_fourier.c
    ../generic_wilson/staggered2naive.c
    wprop_info.c
)

add_milc_executable(lattice_to_scidac lattice_to_scidac.c
    EXTRA_OBJECTS g_objects
    EXTRA_SOURCES
    ../generic/nersc_cksum.c
    ../generic/gauge_utilities.c
    ../generic/check_unitarity.c
    gauge_info_dummy.c
    ../generic/d_plaq4.c
    ../generic/io_lat4.c
)

add_milc_executable(make_modulation make_modulation.c
    DEFINES NO_GAUGE_FIELD
    EXTRA_SOURCES
    ../generic/field_utilities.c
    ../generic/io_scidac.c
    ../generic/io_scidac_types.c
)

add_milc_executable(ckpt1_to_v5 ckpt1_to_v5.c
    EXTRA_SOURCES
    ../generic/gauge_utilities.c
    ../generic/nersc_cksum.c
    ../generic/check_unitarity.c
    gauge_info.c
    ../generic/d_plaq4.c
    ../generic/io_lat4.c
)

add_milc_executable(v5_2scidac v5_to_scidac.c
    DEFINES NO_GAUGE_FIELD
    EXTRA_OBJECTS g_objects
    EXTRA_SOURCES gauge_info_dummy.c
)

add_milc_executable(scidac_2v5 scidac_to_v5.c
    DEFINES NO_GAUGE_FIELD
    EXTRA_OBJECTS g_objects
    EXTRA_SOURCES gauge_info_scidac.c
)

add_milc_executable(crc32 crc32.c)

add_milc_executable(cksum cksum.c)

add_milc_executable(extract_gauge_v5 extract_gauge_v5.c)

add_milc_executable(get_qmp_node_number get_qmp_node_number.c)
