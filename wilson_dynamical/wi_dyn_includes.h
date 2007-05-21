/****************** wi_dyn_includes.h ******************************/
/*
*  Include files for Wilson dynamical application
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
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
Real action();
int congrad_w(int niter,Real rsqmin,Real *final_rsq_ptr);
double d_action();
void gauge_field_copy(field_offset src,field_offset dest);
int congrad(int niter,Real rsqmin,Real *final_rsq_ptr);
int f_measure2();
void grsource();
void checkmul();
void reunitarize();
int s_props( );
int t_props( );
int update();
int w_spectrum();
void boundary_flip(int sign );
int w_spectrum_z();

void init_analyze();
void analyze(int meascount);
void end_analyze(int meascount);

