add_interface_library(my_objects
    gauge_info.c
    ks_source_info.c
    ksprop_info.c
    make_prop.c
    setup.c
    spectrum_ks.c
)


add_interface_library(g_u1_objects
    ../generic_u1/io_u1lat.c
    ../generic_u1/u1link.c
    ../generic_u1/u1plaq.c
)

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
    gauge_objects
    ks_io_objects
    layout
    fft_objects
    fermion_objects
    ks_spectrum_objects
    standard_objects
    eigen_objects
    add_objects
)

add_interface_library(eigcg_objects
    ../generic_ks/inc_eigcg.c
)

add_interface_library(mu_fast_objects
    ../generic_ks/mu_fast.c
)

add_interface_library(naive_debug_objects
    ../generic_wilson/dslash_w3.c
    ../generic_wilson/staggered2naive.c
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
    crypto
)

set(hisq_options
    FN
    HAVE_KS
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
)

set(asqtad_options
    FN
    HAVE_KS
)

set(hisq_gb_baryon_base
    EIGEN_QIO
    GB_BARYON
)

add_milc_executable(ks_spectrum_asqtad control.c
    QUARK asqtad_action
    DEFINES ${asqtad_options}
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(ks_spectrum_eigcg_asqtad control.c
    QUARK asqtad_action
    DEFINES ${asqtad_options} EIGMODE=EIGCG
    EXTRA_OBJECTS fn_objects eigcg_objects
)

add_milc_executable(ks_spectrum_asqtad_eos control.c
    QUARK asqtad_action
    DEFINES ${asqtad_options} DM_DU0 CHEM_POT
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(ks_spectrum_asqtad_mu_eos control.c
    QUARK asqtad_action
    DEFINES ${asqtad_options} TR_MM_INV DM_DU0 CHEM_POT D_CHEM_POT
    EXTRA_OBJECTS fn_objects mu_fast_objects
)

add_milc_executable(ks_spectrum_naive control.c
    QUARK onelink_action
    DEFINES FN HAVE_KS
    EXTRA_OBJECTS fn_objects
)
target_compile_options(ks_spectrum_naive PRIVATE -UD_FN_GATHER13)

add_milc_executable(ks_spectrum_naive_debug control.c
    QUARK onelink_action
    DEFINES FN HAVE_KS DEBUG_NAIVE
    EXTRA_OBJECTS fn_objects naive_debug_objects
)
target_compile_options(ks_spectrum_naive_debug PRIVATE -UD_FN_GATHER13)

add_milc_executable(ks_spectrum_asqtad_u1 control.c
    QUARK asqtad_action
    DEFINES ${asqtad_options} U1_FIELD
    EXTRA_OBJECTS fn_objects g_u1_objects
)

add_milc_executable(ks_spectrum_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options}
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(ks_spectrum_eigcg_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} EIGMODE=EIGCG
    EXTRA_OBJECTS hisq_objects eigcg_objects
)

add_milc_executable(ks_spectrum_hisq_eos control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} DM_DEPS CHEM_POT
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(ks_spectrum_hisq_mu_eos control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} DM_DEPS CHEM_POT D_CHEM_POT TR_MM_INV
    EXTRA_OBJECTS hisq_objects mu_fast_objects
)

add_milc_executable(ks_spectrum_hisq_u1 control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} U1_FIELD
    EXTRA_OBJECTS hisq_objects g_u1_objects
)

add_milc_executable(ks_spectrum_hisq_gb_baryon_no_sink_links control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} GB_BARYON_MMAP NO_SINK_LINKS
    EXTRA_OBJECTS hisq_objects gb_baryon_objects
)

add_milc_executable(ks_spectrum_hisq_gb_baryon_blind_no_sink_links control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} GB_BARYON_MMAP BLIND NO_SINK_LINKS
    EXTRA_OBJECTS hisq_objects gb_baryon_objects
)

add_milc_executable(ks_spectrum_hisq_u1_gb_baryon_blind_no_sink_links control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} GB_BARYON_MMAP BLIND NO_SINK_LINKS U1_FIELD
    EXTRA_OBJECTS hisq_objects gb_baryon_objects g_u1_objects
)

add_milc_executable(ks_spectrum_hisq_gb_baryon_blind_no_sink_links_no_gauss_links control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} GB_BARYON_MMAP BLIND NO_SINK_LINKS NO_GAUSS_SMEAR_LINKS
    EXTRA_OBJECTS hisq_objects gb_baryon_objects
)

add_milc_executable(ks_spectrum_hisq_gb_baryon_blind_one_sided_shift control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} GB_BARYON_MMAP BLIND NO_SINK_LINKS ONE_SIDED_SHIFT_GB
    EXTRA_OBJECTS hisq_objects gb_baryon_objects
)

add_milc_executable(ks_spectrum_hisq_gb_baryon_blind_no_mmap control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${hisq_gb_baryon_base} BLIND NO_SINK_LINKS
    EXTRA_OBJECTS hisq_objects gb_baryon_objects
)
