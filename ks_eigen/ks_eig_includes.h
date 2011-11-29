/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/io_ksprop.h"  /* We really should get rid of this */

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
void measure_chirality(su3_vector *src, double *chirality, int parity);

EXTERN gauge_file *startlat_p;
