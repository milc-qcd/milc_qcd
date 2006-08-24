/****************** cl_dyn_includes.h ******************************/
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
#include "../include/generic_clover.h"
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
int update();
void gauge_action(double *result);
double d_action();
void boundary_flip(int sign );
int f_measure_cl();
int w_spectrum_cl();
int t_props_cl( );
int s_props_cl( );
void make_loop_table2();
void path(int *dir,int *sign,int length);
void single_action(int dir,Real *coeff);
void udadu_mu_nu( field_offset lsrc, field_offset rsrc, field_offset mat, 
		 int mu, int nu, int parity );
void udadu_mat_mu_nu( field_offset matsrc, field_offset matdest, 
		     int mu, int nu );
void grsource_w();
void gauge_field_copy(field_offset src,field_offset dest);
void update_u(Real eps);
void update_h(Real eps);
void mult_sigma_mu_nu( wilson_vector *src, wilson_vector *dest, 
		      int mu, int nu );
int bicongrad_cl(int niter,Real rsqmin,Real *final_rsq_ptr);

/* cl_dyn_includes.h */
