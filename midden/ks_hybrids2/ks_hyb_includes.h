/****************** ks_hyb_includes.h ******************************/
/*
*  Include files for Kogut-Susskind hybrid spectrum application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"

/* prototypes for functions in high level code */

int readin(int prompt);
int  setup();
void update_h( Real eps );
void update_u( Real eps );
void gauge_force( Real eps );
void fermion_force( Real eps );
double gauge_action( );
double hmom_action( );
double fermion_action( );

void smear_links( field_offset src, field_offset dest, Real staple_weight );

complex ploop( void );

int spectrum_hybrids();
void make_field_strength();
int mat_invert_cg2(field_offset src, field_offset dest );
int mat_invert_uml2(field_offset src, field_offset dest );
