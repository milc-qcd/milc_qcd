/****************** wi_hyb_includes.h ******************************/
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
#include "../include/generic_wilson.h"

void boundary_flip( int sign );

int congrad(int niter,Real rsqmin,Real *final_rsq_ptr);
void dslash(field_offset src,field_offset dest,int isign,int parity);
void dslash_special(field_offset src,field_offset dest,int isign,
		    int parity,msg_tag **tag,int is_started);
int f_measure2();
void make_field_strength();
void reunitarize();
int setup();
int readin(int prompt);
void smear_links( field_offset src, field_offset dest );
void smear_links( field_offset src, field_offset dest );
int spectrum_hybrids();
int spectrum_pwave();
int mat_invert_cgilu( field_offset src, field_offset dest );
int mat_invert_cg( field_offset src, field_offset dest );
int mat_invert_bicgilu( field_offset src, field_offset dest );

