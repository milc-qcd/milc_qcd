/****************** gluon_prop_includes.h ******************************/
/*
*  Include files for gluon and (staggered) quark propagator application.
*/

/* Include files */
#include "../include/config.h"	/* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

/* prototypes for functions in high level code */
int  setup();
int readin(int prompt);
void gluon_prop();
int quark_prop();
int quark_renorm();

void gaugefixfft(int gauge_dir, Real accel_param, int max_gauge_iter,
		 Real gauge_fix_tol);
void gaugefixfft_combo(int gauge_dir, Real accel_param, int max_gauge_iter,
		       Real gauge_fix_tol, int nvector,
		       field_offset vector_offset[], int vector_parity[],
		       int nantiherm, field_offset antiherm_offset[], 
		       int antiherm_parity[] );

void dslash_eo( field_offset src, field_offset dest, int parity );
void dslash_eo_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void rephase( int flag );

