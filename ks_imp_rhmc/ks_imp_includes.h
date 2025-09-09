/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

#ifndef KS_IMP_INCLUDES_H_
#define KS_IMP_INCLUDES_H_

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/params_rhmc.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/fermion_links.h"

#ifdef FN
#define dslash_site dslash_fn_site
#define dslash_field dslash_fn_field
#endif
#ifdef EO
#define dslash_site dslash_eo_site
#define dslash_field dslash_eo_field
#endif

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Aggregate time to %s %e\n",(string),dtime);
#else
#define STARTTIME
#define ENDTIME(string)
#endif

/* prototypes for functions in this directory */

/* eo_fermion_force_rhmc.c */
void eo_fermion_force_rhmc( Real eps, params_ratfunc *rf, 
			    su3_vector **multi_x, field_offset phi_off,
			    Real my_rsqmin, int niter, int cg_prec,
			    int ff_prec, fermion_links_t *fl );

/* setup.c */
int setup(void);
int readin(int prompt);

/* update_rhmc.c */

enum int_alg_t { INT_LEAPFROG, INT_OMELYAN, INT_2EPS_3TO1, INT_2EPS_2TO1, 
                 INT_2G1F, INT_3G1F, INT_5G1F, INT_6G1F, INT_4MN4FP, INT_4MN5FV, INT_FOURSTEP, 
		 INT_PLAY };

/* Set default integration algorithm */
#ifndef INT_ALG
#define INT_ALG INT_OMELYAN
#endif

int update(void);
const char *ks_int_alg_opt_chr( void );

/* update_h_rhmc.c */

int update_h_rhmc( Real eps, su3_vector **multi_x );
void update_h_gauge( Real eps );
int update_h_fermion( Real eps, su3_vector **multi_x );
void update_u( Real eps );

/* d_action_rhmc.c */
double d_action_rhmc(su3_vector **multi_x, su3_vector *sumvec );
void gauge_field_copy(field_offset src,field_offset dest);
double fermion_action( su3_vector **multi_x, su3_vector *sumvec );
double hmom_action(void);
void plaquette_action(double *ss_plaq, double *st_plaq);

#endif /* KS_IMP_INCLUDES_H_ */

