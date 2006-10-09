/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
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
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#ifdef HAVE_QDP
#include <qdp.h>
#endif

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
int update();
void make_path_table();
void update_h( Real eps );
void update_h_rhmc( int algorithm_flag, Real eps, su3_vector **multi_x );
void update_h_gauge( int algorithm_flag, Real eps );
void update_h_fermion( int algorithm_flag, Real eps, su3_vector **multi_x );
void update_u( Real eps );
void gauge_force( Real eps );
double d_action_rhmc(su3_vector **multi_x, su3_vector *sumvec );
double imp_gauge_action( );
double hmom_action( );
double fermion_action( su3_vector **multi_x, su3_vector *sumvec );
void ranmom();

// RHMC algorithm stuff
void grsource_imp_rhmc( field_offset dest, Real mass, Real *residues, Real *roots,
  int order, int parity, su3_vector **mult_x, su3_vector *sumvec );
void eo_fermion_force_rhmc( int alg_flag, Real eps, int order, Real mass, Real *residues,
  Real *roots, su3_vector **multi_x, field_offset phi_off );
void eo_fermion_force_rhmc_reverse( int alg_flag, Real eps, int order, Real mass, Real *residues,
  Real *roots, su3_vector **multi_x, field_offset phi_off );
void eo_fermion_force_rhmc_oneterm( Real eps, Real residue, field_offset x_off );
void eo_fermion_force_rhmc_twoterms( Real eps, Real residue1, Real residue2,
	field_offset x1_off, field_offset x2_off );
void eo_fermion_force_rhmc_int( Real eps, Real *residues, su3_vector **xxx, int nterms );

int ks_ratinv(	/* Return value is number of iterations taken */
    field_offset src,   /* source vector (type su3_vector) */
    su3_vector **psim,  /* solution vectors */
    Real mass,          /* quark mass */
    Real *roots,        /* the roots */
    int order,          /* order of rational function approx */
    int niter,          /* maximal number of CG interations */
    Real rsqmin,        /* desired residue squared */
    int parity,         /* parity to be worked on */
    Real *final_rsq_ptr /* final residue squared */
    );
int ks_rateval(
    su3_vector *dest,   /* answer vector */
    field_offset src,   /* source vector (for a_0 term) */
    su3_vector **psim,  /* solution vectors  from multiroot CG */
    Real *residues,     /* the residues */
    int order,          /* order of approximation */
    int parity          /* parity to be worked on */
    );


void hvy_pot( field_offset links );
void f_measure( field_offset phi_off, field_offset xxx_off, Real mass );
void g_measure( void );
void gauge_field_copy(field_offset src,field_offset dest);
void clear_latvec(field_offset v,int parity);
void copy_latvec(field_offset src,field_offset dest,int parity);
void scalar_mult_add_latvec(field_offset src1,field_offset src2,
			     Real scalar,field_offset dest,int parity);
void scalar2_mult_add_su3_vector(su3_vector *a, Real s1, su3_vector *b, 
				 Real s2, su3_vector *c);
void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset dest,int parity);
void scalar_mult_latvec(field_offset src,Real scalar,
			field_offset dest,int parity);

void dslash_eo( field_offset src, field_offset dest, int parity );
void dslash_eo_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void checkmul_imp( field_offset src, Real mass );

void rephase( int flag );
void sym_shift(int dir, field_offset src,field_offset dest);
void zeta_shift(int n, int *d, field_offset src, field_offset dest );

