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

#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
double imp_gauge_action( );
double fermion_action( );
void ranmom();


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

void cleanup_gathers(msg_tag *t1[16],msg_tag *t2[16]);
void cleanup_dslash_temps() ;

void rephase( int flag );
EXTERN gauge_file *startlat_p;
