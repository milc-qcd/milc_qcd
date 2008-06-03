/****************** cl_hyb_includes.h ******************************/
/*
 *  Include files for clover hybrid spectrum application
 */

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/dirs.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/dirs.h"


/**
 **  Function prototypes
 **/

void boundary_flip( int sign );

int f_measure2();
void make_field_strength();

void plaquette(Real *ss_plaq,Real *st_plaq);
void reunitarize();
int setup();
int readin(int prompt);
void smear_links( field_offset src, field_offset dest );
void smear_links( field_offset src, field_offset dest );
int spectrum_hybrids();
int spectrum_pwave();

int initial_set() ;

int mat_invert( field_offset src, field_offset dest );
void copy_site_wilson_vector(field_offset src, field_offset dest);

/**
 **        end of the function prototypes
 **/
