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
    ../generic_wilson/gammas.c
    ../generic_ks/spin_taste_ops.c
    ../generic_ks/shift_field.c
    OBJECTS
    machine_dep_io
    commpkg
)

add_interface_library(objects
    OBJECTS
    my_objects
    gauge_objects
    layout
    fft_objects
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

add_interface_library(eigcg_objects
    ../generic_ks/inc_eigcg.c
)

add_interface_library(deflate_objects)
set(adddefines)
if(HAVE_PRIMME)
    append_interface_library(deflate_objects
        ../generic_ks/eigen_stuff_PRIMME.c
    )
elseif(HAVE_ARPACK)
    append_interface_library(deflate_objects
        ../generic_ks/eigen_stuff_ARPACK.c
    )
else()
    append_interface_library(deflate_objects
        ../generic_ks/eigen_stuff_Ritz.c
    )
    list(APPEND adddefines Kalkreuter_Ritz)
endif()




add_milc_executable(ks_measure_asqtad control.c
    QUARK asqtad_action
    DEFINES FN TR_MM_INV ${adddefines}
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(ks_measure_asqtad_eos control.c
    QUARK asqtad_action
    DEFINES FN DM_DU0 TR_MM_INV CHEM_POT ${adddefines}
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(ks_measure_asqtad_susc control.c
    QUARK asqtad_action
    DEFINES FN CHEM_POT TR_MM_INV ${adddefines}
    EXTRA_OBJECTS fn_objects
)

add_milc_executable(ks_measure_asqtad_mu_eos control.c
    QUARK asqtad_action
    DEFINES FN DM_DU0 TR_MM_INV CHEM_POT D_CHEM_POT ${adddefines}
    EXTRA_OBJECTS fn_objects
    EXTRA_SOURCES ../generic_ks/mu_fast.c
)

set(hisq_options
    FN
    HAVE_KS
    HISQ_REUNIT_ALLOW_SVD
    HISQ_REUNIT_SVD_REL_ERROR=1e-8
    HISQ_REUNIT_SVD_ABS_ERROR=1e-8
)

add_milc_executable(ks_measure_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} ${adddefines}
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(ks_measure_eigcg_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} EIGMODE=EIGCG ${adddefines}
    EXTRA_OBJECTS hisq_objects eigcg_objects
)

add_milc_executable(ks_measure_current_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} CURRENT_DISC ${adddefines}
    EXTRA_OBJECTS hisq_objects
    EXTRA_SOURCES ../generic_ks/f_meas_current.c
)

add_milc_executable(ks_measure_current_eigcg_hisq control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} CURRENT_DISC EIGMODE=EIGCG ${adddefines}
    EXTRA_OBJECTS hisq_objects eigcg_objects
    EXTRA_SOURCES ../generic_ks/f_meas_current.c
)

add_milc_executable(ks_measure_hisq_eos control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} DM_DEPS CHEM_POT ${adddefines}
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(ks_measure_hisq_susc control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} CHEM_POT ${adddefines}
    EXTRA_OBJECTS hisq_objects
)

add_milc_executable(ks_measure_hisq_mu_eos control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} DM_DEPS CHEM_POT D_CHEM_POT ${adddefines}
    EXTRA_OBJECTS hisq_objects
    EXTRA_SOURCES ../generic_ks/mu_fast.c
)

add_milc_executable(ks_measure_current_hisq_u1 control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} CURRENT_DISC U1_FIELD ${adddefines}
    EXTRA_OBJECTS hisq_objects g_u1_objects
    EXTRA_SOURCES ../generic_ks/f_meas_current.c
)

add_milc_executable(ks_measure_current_hisq_u1_loop control.c
    QUARK hisq/hisq_u3_action
    DEFINES ${hisq_options} CURRENT_DISC U1_FIELD QED_LOOP ${adddefines}
    EXTRA_OBJECTS hisq_objects g_u1_objects
    EXTRA_SOURCES ../generic_ks/f_meas_current.c
)
