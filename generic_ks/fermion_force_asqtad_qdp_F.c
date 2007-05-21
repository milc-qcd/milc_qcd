/***************** fermion_force_asqtad_qdp_F.c ****************************/
/* MILC Version 7 */

/* For specifically single precision calculations using QDP routines */

/* 5/09/07 C. DeTar */

#undef QDP_Precision
#define QDP_Precision 'F'

#include "fermion_force_asqtad_qdp_P.c"

