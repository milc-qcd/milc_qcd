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
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "params_rhmc.h"
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
		 INT_2G1F, INT_3G1F, INT_4MN4FP, INT_4MN5FV, INT_FOURSTEP, 
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

/* grsource_rhmc.c */

void grsource_imp_rhmc( field_offset dest, params_ratfunc *rf,
			int parity, su3_vector **multi_x, su3_vector *sumvec,
			Real my_rsqmin, int my_niter, int my_prec,
			imp_ferm_links_t *fn, int naik_term_epsilon_index,
			Real naik_term_epsilon);

/* ks_ratinv.c */

int ks_ratinv(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *roots,	/* the roots */
    int order,		/* order of rational function approx */
    int my_niter,	/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precicion */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    imp_ferm_links_t *fn,     /* Fermion links */
    int naik_term_epsilon_index, /* Index of naik term common to this set */
    Real naik_term_epsilon /* Epsilon common to this set */
		);

int ks_rateval(
    su3_vector *dest,   /* answer vector */
    field_offset src,   /* source vector (for a_0 term) */
    su3_vector **psim,  /* solution vectors  from multiroot CG */
    Real *residues,     /* the residues */
    int order,          /* order of approximation */
    int parity          /* parity to be worked on */
    );

/* load_rhmc_params */

params_rhmc *load_rhmc_params(char filename[], int n_pseudo);

/* d_action_rhmc.c */
double d_action_rhmc(su3_vector **multi_x, su3_vector *sumvec );
void gauge_field_copy(field_offset src,field_offset dest);
double fermion_action( su3_vector **multi_x, su3_vector *sumvec );
double hmom_action(void);

#endif /* KS_IMP_INCLUDES_H_ */

