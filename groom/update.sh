#! /bin/sh

# Change to "normal_exit"


FILES="\
 arb_dirac_eigen/setup_p_cl.c \
 arb_dirac_invert/setup_p_cl.c \
 clover_dynamical/setup.c \
 clover_hybrids/setup.c \
 clover_invert/setup_cl.c \
 dense_static_su3/setup.c \
 generic/com_mpi.c \
 generic/io_romio.c \
 h_dibaryon/setup_H_cl.c \
 heavy/setup_mr.c \
 hqet_heavy_to_light/setup_hqet_form.c \
 hvy_qpot/setup.c \
 ks_dynamical/setup.c \
 ks_eigen/setup.c \
 ks_hybrids2/setup.c \
 ks_imp_dyn1/setup.c \
 ks_imp_dyn2/setup.c \
 propagating_form_factor/setup_form.c \
 pure_gauge/setup.c \
 schroed_cl_inv/setup_cl.c \
 schroed_ks_dyn/setup.c \
 schroed_pg/setup.c \
 smooth_inst/setup.c \
 string_break/setup.c \
 symanzik_sl32/setup.c \
 wilson_dynamical/setup.c \
 wilson_hybrids/setup.c \
 wilson_invert/setup_w.c \
 wilson_static/setup_mr.c"

FILES=arb_dirac_eigen/setup_p_cl.c

for file in $FILES
do
  mv $file $file.bak
sed 's/  if( par_buf.stopflag != 0 ){\
\#ifdef MPI\
    MPI_Finalize();\
\#endif\
    exit(0);\
  }/  if( par_buf.stopflag != 0 )\
    normal_exit(0)/' $file.bak > $file
done