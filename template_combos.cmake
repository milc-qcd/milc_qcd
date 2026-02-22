# This is not a standalone Makefile
# It is intended to be included in other Makefiles

# It defines standard ccmbinations of generic object files for building code
# These definitions are intended to apply to all builds.
# If you want to customize a build, this list can be used as a starting point
# for making changes.

######################################################################
# Standard lattice routines
######################################################################

# These are required for essentially all applications

add_interface_library(standard_objects
    generic/blind_data.c
    generic/field_translation.c
    generic/field_utilities.c
    generic/gauge_utilities.c
    generic/io_detect.c
    generic/io_helpers.c
    generic/io_lat_utils.c
    generic/make_lattice.c
    generic/mmap_cache.c
    generic/ranstuff.c
    generic/remap_stdio_from_args.c
)
if(HAVEQIO)
    append_interface_library(standard_objects
        generic/file_types_milc_usqcd.c generic/io_scidac.c generic/io_scidac_types.c
    )
endif()
if(HAVE_APE_IO)
    append_interface_library(standard_objects
        generic/io_ape_links.c
    )
endif()

# These use the gauge field extensively

add_interface_library(gauge_objects
    generic/ape_smear.c
    generic/check_unitarity.c
    generic/d_plaq4.c
    generic/gaugefix2.c
    generic/io_lat4.c
    generic/momentum_twist.c
    generic/nersc_cksum.c
    generic/path_product.c
    generic/project_su3_hit.c
    generic/reunitarize2.c
    generic/show_generic_opts.c
    generic/show_scidac_opts.c
)
if(HAVE_GRID)
    append_interface_library(gauge_objects
        generic/gridMap.cc generic/milc_to_grid_utilities.cc
    )
endif()
if(HAVEQOP)
    append_interface_library(gauge_objects
        generic/map_milc_to_qopqdp.c generic/milc_to_qop_utilities.c
    )
endif()
if(HAVE_QUDA)
    append_interface_library(gauge_objects
        generic/milc_to_quda_utilities.c generic/d_plaq4_gpu.c generic/ploop3_gpu.c
    )
endif()
if(HAVE_QPHIX)
    append_interface_library(gauge_objects
        generic/map_milc_to_qphix.c generic/milc_to_qphix_utilities.c
    )
endif()
if(HAVE_QPHIXJ)
    append_interface_library(gauge_objects
        generic/qphixjClovMap.cc generic/milc_to_qphixj_utilities.cc
    )
endif()

if(HAVEFFTW)
    add_interface_library(fft_objects
        generic/remap_fftw_fourier.c
    )
else()
    add_interface_library(fft_objects
        generic/restrict_fourier.c
    )
endif()

######################################################################
# Standard fermion routines
######################################################################

# These are needed for either KS or Clover sources:

add_interface_library(fermion_objects
    generic/discretize_wf.c
    generic_wilson/gammas.c
    generic/io_source_cmplx_fm.c
    generic/phases.c
    generic/quark_source.c
    generic/quark_source_io.c
    generic/quark_source_sink_op.c
    generic_ks/spin_taste_ops.c
    generic_ks/shift_field.c
)

# These are used to build propagators

add_interface_library(ks_links_objects
    generic_ks/charge_utilities.c
    generic_ks/fermion_links_from_site.c
    generic_ks/f_meas.c
    generic_ks/gauss_smear_ks.c
    generic_ks/gauss_smear_ks_cpu.c
    generic_ks/gauss_smear_ks_QUDA.c
    generic_ks/grsource_imp.c
    generic_ks/naik_eps_utilities.c
    generic_ks/path_transport.c
    generic_ks/rephase.c
    generic_ks/show_generic_ks_opts.c
    generic_ks/show_hisq_links_opts.c
)

# These are used to read propagators

add_interface_library(ks_io_objects
    generic_ks/io_helpers_ks.c
    generic_ks/io_prop_ks.c
    generic_ks/io_prop_ks_fm.c
)
if(HAVEQIO)
    append_interface_library(ks_io_objects
        generic_ks/io_scidac_ks.c
    )
endif()

# These are used to compute the hadron correlators

add_interface_library(ks_spectrum_objects
    generic_ks/ks_baryon.c
)
if(HAVE_KS_CONT_GPU)
    # use QUDA GPU contraction code
    append_interface_library(ks_spectrum_objects
        generic_ks/ks_meson_mom_quda.c
    )
else()
    # otherwise use old CPU contraction code
    append_interface_library(ks_spectrum_objects
        generic_ks/ks_meson_mom.c
    )
endif()

add_interface_library(gb_baryon_objects
    generic_ks/gb_baryon_mmap.c
    generic_ks/gb_baryon_snk.c
    generic_ks/gb_baryon_src.c
    generic_ks/gb_baryon_3pt.c
    generic_ks/gb_ops.c
)

# These are used to build propagators

add_interface_library(cl_links_objects
    generic_wilson/gauss_smear_w.c
)

# These are used to read propagators

add_interface_library(cl_io_objects
    generic_wilson/canopy2weyl_rot.c
    generic_wilson/io_helpers_w.c
    generic_wilson/io_prop_w.c
    generic_wilson/staggered2naive.c
)
if(HAVEQIO)
    append_interface_library(cl_io_objects
        generic_wilson/io_scidac_w.c
    )
endif()

# These are used to compute the spectrum

add_interface_library(cl_spectrum_objects
    generic_wilson/baryon_cont.c
    generic_wilson/w_baryon.c
    generic_wilson/w_baryon_hl.c
    generic_wilson/w_meson_mom.c
    generic_wilson/w_meson_open_mom.c
)

# These are used to calculate eigenvectors

add_interface_library(eigen_objects
    generic_ks/eigen_stuff_helpers.c generic_ks/io_helpers_ks_eigen.c generic_ks/io_ks_eigen.c generic_ks/read_eigen_param.c
)
if(HAVEQIO)
    append_interface_library(eigen_objects
        generic_ks/io_scidac_ks_eigen.c generic_ks/io_grid_ks_eigen.c
    )
endif()
if(HAVE_EIG_GPU)
    if(HAVE_GRID)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_Grid.c generic_ks/gridStaggEigen.cc
        )
    elseif(HAVE_QUDA)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_QUDA.c
        )
    else()
        append_interface_library(eigen_objects
            generic_ks/must_specify_HAVE_QUDA_or_HAVE_GRID.c
        )
    endif()
else()
    if(HAVE_PRIMME)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_PRIMME.c
        )
    elseif(HAVE_QDP)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_qdp.c
        )
    elseif(HAVE_ARPACK)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_ARPACK.c
        )
    elseif(HAVEQDP)
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_qdp.c generic_ks/jacobi.c
        )
    else()
        # Default solver
        append_interface_library(eigen_objects
            generic_ks/eigen_stuff_Ritz.c generic_ks/jacobi.c
        )
    endif()
endif()

######################################################################
# Staggered fermion links routines
######################################################################

# Standard MILC combinations

add_interface_library(flinks
    generic_ks/fermion_links.c generic_ks/fermion_links_fn_load_milc.c
    generic_ks/fermion_links_fn_twist_milc.c generic/general_staple.c generic_ks/fn_links_milc.c
)

add_interface_library(flinks_fn_milc
    generic_ks/fermion_links_milc.c generic_ks/ks_action_paths.c
    OBJECTS flinks
)

add_interface_library(flinks_hisq_milc
    generic_ks/fermion_links_hisq_milc.c generic_ks/fermion_links_hisq_load.c
    generic_ks/ks_action_paths_hisq.c generic_ks/su3_mat_op.c generic/stout_smear.c
    OBJECTS flinks
)

add_interface_library(flinks_eo_milc
    generic_ks/fermion_links.c generic_ks/fermion_links_milc.c generic_ks/fermion_links_eo_load_milc.c
    generic_ks/eo_links_milc.c generic_ks/ks_action_paths.c
)

# QUDA support

add_interface_library(flinks_fn_quda
    generic_ks/fermion_links_milc.c generic_ks/fermion_links_fn_load_quda.c
    generic_ks/ks_action_paths.c
    OBJECTS flinks
)

add_interface_library(flinks_hisq_quda
    generic_ks/fermion_links.c generic_ks/fermion_links_fn_load_quda.c
    generic_ks/fermion_links_fn_twist_milc.c generic_ks/fn_links_milc.c
    generic_ks/ks_action_paths_hisq.c generic_ks/su3_mat_op.c generic/stout_smear.c
    generic_ks/fermion_links_hisq_milc.c generic_ks/fermion_links_hisq_load.c
)

# GRID support

add_interface_library(flinks_fn_grid
    generic_ks/fermion_links_milc.c
    generic_ks/fermion_links_hisq_load_grid.c generic_ks/fermion_links_hisq_load_grid_D.c generic_ks/fermion_links_hisq_load_grid_F.c
    generic_ks/gridHISQLinks.cc generic_ks/ks_action_paths.c
    OBJECTS flinks
)

add_interface_library(flinks_hisq_grid
    generic_ks/fermion_links_hisq_milc.c
    generic_ks/fermion_links_hisq_load.c
    generic_ks/fermion_links_hisq_load_grid.c generic_ks/fermion_links_hisq_load_grid_D.c generic_ks/fermion_links_hisq_load_grid_F.c
    generic_ks/gridHISQLinks.cc generic_ks/ks_action_paths_hisq.c generic_ks/su3_mat_op.c generic/stout_smear.c
    OBJECTS flinks
)

# Standard QOP combinations

add_interface_library(flinks_qop
    generic_ks/fermion_links.c
)

add_interface_library(flinks_fn_qop
    generic_ks/fn_links_qop.c generic_ks/fermion_links_asqtad_qop.c
    generic_ks/fermion_links_fn_twist_qop.c generic_ks/ks_action_coeffs_asqtad_qop.c generic_ks/ks_action_paths.c
    OBJECTS flinks_qop
)

add_interface_library(flinks_hisq_qop
    generic_ks/fn_links_qop.c generic_ks/hisq_links_qop.c
    generic_ks/fermion_links_hisq_qop.c generic_ks/fermion_links_fn_twist_qop.c
    generic_ks/ks_action_coeffs_hisq_qop.c generic_ks/ks_action_paths_hisq.c
    generic_ks/su3_mat_op.c
    OBJECTS flinks_qop
)

# Generic actions are not supported in QOP, so we use MILC
add_interface_library(flinks_eo_qop OBJECTS flinks_eo_milc)

if(HAVE_FL_GPU)
    if(HAVE_QUDA)
        add_interface_library(flinks_fn OBJECTS flinks_fn_quda)
        add_interface_library(flinks_eo OBJECTS flinks_eo_milc)
        add_interface_library(flinks_hisq OBJECTS flinks_hisq_quda)
    elseif(HAVE_GRID)
        add_interface_library(flinks_fn OBJECTS flinks_fn_grid)
        add_interface_library(flinks_eo OBJECTS flinks_eo_milc)
        add_interface_library(flinks_hisq OBJECTS flinks_hisq_grid)
    endif()
else()
    if(HAVEQOP)
        add_interface_library(flinks_fn OBJECTS flinks_fn_qop)
        add_interface_library(flinks_eo OBJECTS flinks_eo_qop)
        add_interface_library(flinks_hisq OBJECTS flinks_hisq_qop)
    else()
        add_interface_library(flinks_fn OBJECTS flinks_fn_milc)
        add_interface_library(flinks_eo OBJECTS flinks_eo_milc)
        add_interface_library(flinks_hisq OBJECTS flinks_hisq_milc)
    endif()
endif()

######################################################################
# Staggered Dslash routine
######################################################################

# Standard MILC

# Choices here are dslash_fn.o dslash_fn2.o dslash_fn_dblstore.o
if(HAVE_QUDA)
    # When using QUDA, the back links are not used and just add unnecessary overhead
    add_interface_library(dslash_fn_milc
        generic_ks/dslash_fn.c
    )
else()
    add_interface_library(dslash_fn_milc
        generic_ks/dslash_fn_dblstore.c
    )
endif()

# No other choice
add_interface_library(dslash_eo
    generic_ks/dslash_eo.c
)

add_interface_library(dslash_cl
    generic_wilson/dslash_w_space.c generic_wilson/dslash_w3.c
)

# Standard QOP

add_interface_library(dslash_fn_qop
    generic_ks/dslash_fn_qop.c
)

######################################################################
# Staggered single-mass inverters
######################################################################

# Standard MILC combinations

add_interface_library(congrad_fn_base
    generic_ks/mat_invert.c generic_ks/ks_invert.c generic_ks/d_congrad5_fn.c generic_ks/d_congrad_opt.c generic/report_invert_status.c
)

add_interface_library(congrad_fn_milc_cpu
    generic_ks/d_congrad5_two_src.c generic_ks/d_congrad5_fn_milc.c
    OBJECTS congrad_fn_base
)

# No other choice
add_interface_library(congrad_eo
    generic_ks/d_congrad5_eo.c generic_ks/d_congrad_opt.c generic_ks/mat_invert.c generic_ks/ks_invert.c generic/report_invert_status.c
)

# GRID support

add_interface_library(congrad_fn_grid
    generic_ks/d_congrad5_two_src.c generic_ks/d_congrad5_fn_grid.c
    generic_ks/d_congrad5_fn_grid_D.c generic_ks/d_congrad5_fn_grid_F.c generic_ks/d_congrad5_fn_milc.c
    generic_ks/gridStaggInvert.cc
    OBJECTS congrad_fn_base
)

# QUDA support

add_interface_library(congrad_fn_quda
    generic_ks/d_congrad5_fn_quda.c
    OBJECTS congrad_fn_milc_cpu
)

# QPHIX support

add_interface_library(congrad_fn_qphix
    generic_ks/d_congrad5_two_src.c generic_ks/d_congrad5_fn_qphix.c
    generic_ks/d_congrad5_fn_qphix_D.c generic_ks/d_congrad5_fn_qphix_F.c generic_ks/d_congrad5_fn_milc.c
    OBJECTS congrad_fn_base
)

# Standard QOP combinations

add_interface_library(congrad_fn_qop
    generic_ks/d_congrad5_fn_qop_two_src.c generic_ks/d_congrad5_fn_qop.c
    generic_ks/d_congrad5_fn_qop_D.c generic_ks/d_congrad5_fn_qop_F.c
    OBJECTS congrad_fn_base
)

if(HAVE_FN_CG_GPU)
    if(HAVE_QUDA)
        add_interface_library(congrad_fn_milc OBJECTS congrad_fn_quda)
    elseif(HAVE_GRID)
        add_interface_library(congrad_fn_milc OBJECTS congrad_fn_grid)
    endif()
else()
    if(HAVE_QPHIX)
        add_interface_library(congrad_fn_milc OBJECTS congrad_fn_qphix)
    else()
        add_interface_library(congrad_fn_milc OBJECTS congrad_fn_milc_cpu)
    endif()
endif()

######################################################################
# Staggered multimass inverters
######################################################################

# Standard MILC combinations

add_interface_library(multi_inv_fn_milc_cpu
    generic_ks/ks_multicg.c generic_ks/ks_multicg_offset.c
)

# No other choice
add_interface_library(multi_inv_eo
    generic_ks/ks_multicg.c generic_ks/ks_multicg_offset.c
)

# QUDA support

add_interface_library(multi_inv_fn_quda
    generic_ks/ks_multicg_offset_quda.c
    OBJECTS multi_inv_fn_milc_cpu
)

# Grid support

add_interface_library(multi_inv_fn_grid
    generic_ks/ks_multicg_offset_grid.c
    generic_ks/ks_multicg_offset_grid_D.c generic_ks/ks_multicg_offset_grid_F.c
    OBJECTS multi_inv_fn_milc_cpu
)

# QPHIX support

add_interface_library(multi_inv_fn_qphix
    generic_ks/ks_multicg.c generic_ks/ks_multicg_offset_qphix.c
    generic_ks/ks_multicg_offset_qphix_D.c generic_ks/ks_multicg_offset_qphix_F.c
)

# Standard QOP combinations

add_interface_library(multi_inv_fn_qop
    generic_ks/ks_multicg.c generic_ks/ks_multicg_offset_qop.c
    generic_ks/ks_multicg_offset_qop_D.c generic_ks/ks_multicg_offset_qop_F.c
)

if(HAVE_FN_CG_GPU)
    if(HAVE_QUDA)
        add_interface_library(multi_inv_fn_milc OBJECTS multi_inv_fn_quda)
    elseif(HAVE_GRID)
        add_interface_library(multi_inv_fn_milc OBJECTS multi_inv_fn_grid)
    endif()
else()
    if(HAVE_QPHIX)
        add_interface_library(multi_inv_fn_milc OBJECTS multi_inv_fn_qphix)
    else()
        add_interface_library(multi_inv_fn_milc OBJECTS multi_inv_fn_milc_cpu)
    endif()
endif()

######################################################################
# Staggered fermion force routines
######################################################################

# Standard MILC CPU combinations

add_interface_library(asq_force_milc
    generic_ks/fermion_force_asqtad.c generic_ks/fermion_force_multi.c
    generic_ks/fermion_force_fn_multi.c generic_ks/ff_opt.c
)

add_interface_library(hisq_force_milc
    generic_ks/fermion_force_hisq_multi.c generic_ks/fermion_force_hisq_multi_cpu.c
    generic_ks/show_hisq_force_opts.c
)

# For 2 and 2+1 flavor (one-term and two-term only)

add_interface_library(eo_force
    generic_ks/fermion_force_eo_milc.c generic_ks/show_hisq_force_opts.c
)

# QOP support

add_interface_library(asq_force_qop
    generic_ks/fermion_force_asqtad_qop.c generic_ks/fermion_force_asqtad_qop_F.c
    generic_ks/fermion_force_asqtad_qop_D.c generic_ks/ff_opt.c # ${FORCE_OPTS}
)

add_interface_library(hisq_force_qop
    generic_ks/fermion_force_hisq_qop.c generic_ks/fermion_force_hisq_qop_F.c
    generic_ks/fermion_force_hisq_qop_D.c generic_ks/show_hisq_force_opts.c
)

# QUDA support

add_interface_library(asq_force_quda
    generic_ks/fermion_force_asqtad_gpu.c
    OBJECTS asq_force_milc
)

add_interface_library(hisq_force_quda
    generic_ks/fermion_force_hisq_multi.c generic_ks/fermion_force_hisq_multi_quda.c
    generic_ks/show_hisq_force_opts.c
)

# GRID support

add_interface_library(hisq_force_grid
    generic_ks/fermion_force_hisq_multi.c generic_ks/fermion_force_hisq_multi_grid.c
    generic_ks/show_hisq_force_opts.c
)

# Define ASQ_FORCE and HISQ_FORCE depending on compilation parameters

if(HAVE_FF_GPU)
    if(HAVE_QUDA)
        add_interface_library(asq_force OBJECTS asq_force_quda)
        add_interface_library(hisq_force OBJECTS hisq_force_quda)
    elseif(HAVE_GRID)
        add_interface_library(hisq_force OBJECTS hisq_force_grid)
    endif()
else()
    if(HAVEQOP)
        add_interface_library(asq_force OBJECTS asq_force_qop)
        add_interface_library(hisq_force OBJECTS hisq_force_qop)
    else()
        add_interface_library(asq_force OBJECTS asq_force_milc)
        add_interface_library(hisq_force OBJECTS hisq_force_milc)
    endif()
endif()

######################################################################
# Clover link creation
######################################################################

add_interface_library(flinks_cl
    generic_clover/f_mu_nu.c generic_clover/make_clov2.c
)

######################################################################
# Clover inverters
######################################################################

add_interface_library(congrad_cl_base
    generic_clover/d_cgilu_cl.c generic_clover/d_hopilu_cl.c generic_clover/d_mrilu_cl.c generic_clover/cl_solver_utilities.c
    generic_wilson/wilson_invert.c
)

add_interface_library(congrad_cl_qop
    generic_clover/d_bicgilu_cl_qop.c generic_clover/d_bicgilu_cl_qop_D.c generic_clover/d_bicgilu_cl_qop_F.c
    OBJECTS congrad_cl_base
)

# Standard MILC combinations
add_interface_library(congrad_cl_milc_cpu
    generic_clover/d_bicgilu_cl.c
    OBJECTS congrad_cl_base
)

# GPU support
add_interface_library(congrad_cl_milc_gpu
    generic_clover/d_bicgilu_cl_gpu.c
    OBJECTS congrad_cl_base
)

# QPHIXJ support
# add_interface_library(congrad_cl_qphixj
#     generic_clover/d_bicgilu_cl_qphixj.c generic_clover/d_bicgilu_cl_qphixj_F.c
#     generic_clover/d_bicgilu_cl_qphixj_D.c
#     OBJECTS congrad_cl_base
# )

# QPHIXJ support
add_interface_library(congrad_cl_milc_qphixj
    generic_clover/d_bicgilu_cl_qphixj.c generic_clover/qphixjClovInvert.cc
    generic_clover/d_bicgilu_cl_qphixj_F.c generic_clover/d_bicgilu_cl_qphixj_D.c
    OBJECTS congrad_cl_base
)

if(HAVE_CL_GPU)
    add_interface_library(congrad_cl_milc OBJECTS congrad_cl_milc_gpu)
elseif(HAVE_QPHIXJ)
    add_interface_library(congrad_cl_milc OBJECTS congrad_cl_milc_qphixj)
else()
    add_interface_library(congrad_cl_milc OBJECTS congrad_cl_milc_cpu)
endif()

######################################################################
# Gauge force routines
######################################################################

# Standard MILC combinations

add_interface_library(gauge_force_milc_cpu
    generic/gauge_force_imp.c generic/gauge_action_imp.c generic/gauge_stuff.c generic/ranmom.c
)

# GPU support

add_interface_library(gauge_force_milc_gpu
    generic/gauge_force_imp_gpu.c generic/gauge_action_imp_gpu.c generic/gauge_action_imp.c generic/gauge_stuff.c generic/ranmom.c
)

# Standard QOP combinations

add_interface_library(gauge_force_qop
    generic/gauge_force_symzk1_qop.c generic/gauge_action_imp.c generic/gauge_stuff.c generic/ranmom.c
)

# Standard QPhiX combinations

add_interface_library(gauge_force_qphix
    generic/gauge_force_symzk1_qphix.c generic/gauge_action_imp.c generic/gauge_stuff.c generic/ranmom.c
    generic/gauge_force_symzk1_qphix_D.c generic/gauge_force_symzk1_qphix_F.c
)

if(HAVE_GF_GPU)
    add_interface_library(gauge_force_milc OBJECTS gauge_force_milc_gpu)
else()
    if(HAVE_GF_QPHIX)
        add_interface_library(gauge_force_milc OBJECTS gauge_force_qphix)
    else()
        add_interface_library(gauge_force_milc OBJECTS gauge_force_milc_cpu)
    endif()
endif()

######################################################################
# QOP or MILC/GPU or QPHIX or GRID
######################################################################

if(HAVEQOP)
    # Interface to access QOP
    add_interface_library(congrad_cl OBJECTS congrad_cl_qop)
    add_interface_library(dslash_fn OBJECTS dslash_fn_qop)
    add_interface_library(congrad_fn OBJECTS congrad_fn_qop)
    add_interface_library(multi_inv_fn OBJECTS multi_inv_fn_qop)
    add_interface_library(gauge_force OBJECTS gauge_force_qop)
else()
    add_interface_library(congrad_cl OBJECTS congrad_cl_milc)
    add_interface_library(dslash_fn OBJECTS dslash_fn_milc)
    add_interface_library(congrad_fn OBJECTS congrad_fn_milc)
    add_interface_library(multi_inv_fn OBJECTS multi_inv_fn_milc)
    add_interface_library(gauge_force OBJECTS gauge_force_milc)
endif()

######################################################################
# Standard lists of objects for staggered fermion links and inverters
######################################################################

add_interface_library(cl_objects
    OBJECTS
    cl_links_objects
    flinks_cl
    congrad_cl
    dslash_cl
)

add_interface_library(eo_objects
    OBJECTS
    ks_links_objects
    flinks_eo
    congrad_eo
    dslash_eo
    multi_inv_eo
)

add_interface_library(fn_objects
    OBJECTS
    ks_links_objects
    flinks_fn
    congrad_fn
    dslash_fn
    multi_inv_fn
)

add_interface_library(hisq_objects
    OBJECTS
    ks_links_objects
    flinks_hisq
    congrad_fn
    dslash_fn
    multi_inv_fn
)
