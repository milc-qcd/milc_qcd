/* ****************************************************
   Intended for debugging features
   AB 1/28/8, started
   **************************************************** */

#ifndef KS_IMP_RHMC_DEBUG_H

#include "ks_imp_includes.h"    /* definitions files and prototypes */
//#include "../generic_ks/generic_ks_includes.h"
#include "lattice_qdp.h"

/* calculate plaquettes and write into a text file */
void d_plaquette_field_dump(su3_matrix **U_field, char *file_name_prefix);

/* Measure plaquettes on fat links */
void g_measure_plaq();

/* Measure plaquettes and write into file */
void g_measure_tune();


#define KS_IMP_RHMC_DEBUG_H

#endif /* KS_IMP_RHMC_DEBUG_H */

