add_interface_library(my_objects
    gauge_info.c
    setup.c
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
    ../generic_ks/spin_taste_ops.c
    ../generic_ks/shift_field.c
    ../generic_wilson/gammas.c
    OBJECTS
    machine_dep_io
    commpkg
)

add_interface_library(objects
    OBJECTS
    my_objects
    fft_objects
    gauge_objects
    layout
    standard_objects
    eigen_objects
    add_objects
)

add_interface_library(libraries
    OBJECTS
    qcd
    scidac
    add
    m
)

add_milc_executable(su3_eigen_asqtad control.c
    QUARK asqtad_action
    DEFINES FN EIGEN_QIO
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(su3_eigen_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES FN EIGEN_QIO
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(su3_eigen_hisq_pack control.c
    QUARK hisq/hisq_u3_action
    DEFINES FN EIGEN_QIO PACK_EIGEN
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(su3_u1_eigen_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES FN EIGEN_QIO U1_FIELD
    EXTRA_OBJECTS hisq_objects g_u1_objects
)

add_milc_executable(su3_eigen_p4fat3 control.c
    QUARK p4fat3_action
    DEFINES EO
    EXTRA_OBJECTS eo_objects
)
target_compile_options(su3_eigen_p4fat3 PRIVATE -UKS_MULTICG -UHAVE_QOP -UHAVE_QDP)
