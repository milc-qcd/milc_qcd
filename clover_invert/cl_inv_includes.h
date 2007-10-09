/****************** cl_inv_includes.h ******************************/
/*
*  Include files for the clover_invert application
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
#include "../include/comdefs.h"	/* definitions and variables for communications */
#include "../include/io_ksprop.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int setup_cl();
int readin(int prompt);

/* dirac_utilities.c */
double start_timing(void);
void print_timing(double dtime, char *str);
spin_wilson_vector *create_swv_field(void);
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field wp, int color);
wilson_prop_field create_wp_field(void);
wilson_vector *create_wv_field(void);
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field wp, 
		     int color, int spin);
void copy_wp_from_wv(wilson_prop_field wp, wilson_vector *wv, 
		     int color, int spin);
void destroy_swv_field(spin_wilson_vector *swv);
void destroy_wv_field(wilson_vector *wv);
void destroy_wp_field(wilson_prop_field wp);

/* spectrum_cl.c */
void spectrum_cl_init();
void spectrum_cl_baryon(wilson_prop_field qp, complex *bp[]);
void rotate_prop(spin_wilson_vector *rp, wilson_prop_field qp, int color);
void sink_smear_prop(wilson_prop_field qp);
void spectrum_cl_diag_gen_meson(wilson_prop_field qp, complex *mp[]);
void spectrum_cl_diag_meson(wilson_prop_field qp);
void spectrum_cl_diag_rot_meson(wilson_prop_field qp);
void spectrum_cl_diag_smeared_meson(wilson_prop_field qp);
void spectrum_cl_print(int t0, int k);
void spectrum_cl_cleanup();
void spectrum_cl(wilson_prop_field qp, int t0, int k);

/* spectrum_cl_hl.c */
void spectrum_cl_hl_init();
void spectrum_cl_hl_diag_baryon(wilson_prop_field qp, int k);
void spectrum_cl_hl_diag_meson(wilson_prop_field qp, int k);
void spectrum_cl_hl_diag_rot_meson(wilson_prop_field qp, int k);
void spectrum_cl_hl_diag_smeared_meson(wilson_prop_field qp, int k);
void spectrum_cl_hl_baryon(wilson_prop_field qp1, wilson_prop_field qp2,
			   complex *bpa[], complex *bpb[]);
void spectrum_cl_hl_offdiag_baryon(wilson_prop_field qp1, 
				   wilson_prop_field qp2, 
				   int j, int k);
void spectrum_cl_hl_offdiag_gen_meson(wilson_prop_field qp1, 
				      wilson_prop_field qp2, 
				      complex *mp[]);
void spectrum_cl_hl_offdiag_meson(wilson_prop_field qp1, 
				  wilson_prop_field qp2, 
				  int j, int k);
void spectrum_cl_hl_offdiag_rot_meson(wilson_prop_field qp1, 
				      wilson_prop_field qp2, 
				      int j, int k);
void spectrum_cl_hl_offdiag_smeared_meson(wilson_prop_field qp1, 
					  wilson_prop_field qp2, 
					  int j, int k);
void spectrum_cl_hl_print(int t0);
void spectrum_cl_hl_cleanup();

/*  cl_inv_includes.h */




