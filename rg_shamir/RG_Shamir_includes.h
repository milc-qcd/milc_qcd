/****************** RG_Shamir_includes.h ******************************/
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
void update_u( Real eps );
void gauge_force( Real eps );
double imp_gauge_action( );
double hmom_action( );
double fermion_action( );
void ranmom();


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

