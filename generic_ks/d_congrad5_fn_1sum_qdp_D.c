/******************** d_congrad5_fn_1sum_qdp_D.c ****************************/
/* MILC Version 7 */

/* Dor specifically single precision inversions using QDP routines */

/* 12/03/06 C. DeTar */

#undef QDP_Precision
#define QDP_Precision 'D'

#define KS_CONGRAD_QDP       ks_congrad_qdp_D
#define KS_CONGRAD_MILC2QDP  ks_congrad_milc2qdp_D

#include "d_congrad5_fn_1sum_qdp_P.c"

