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

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Time to %s %e\n",(string),dtime);
#else
#define STARTTIME
#define ENDTIME(string)
#endif

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);
int num_mes_report(void);

/* dirac_utilities.c */
double start_timing(void);
void print_timing(double dtime, char *str);
spin_wilson_vector *create_swv_field(void);
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field wp, int color);
wilson_prop_field create_wp_field(void);
wilson_vector *create_wv_field(void);
su3_vector *create_v_field(void);
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field wp, 
		     int color, int spin);
void copy_wp_from_wv(wilson_prop_field wp, wilson_vector *wv, 
		     int color, int spin);
void copy_wp_field(wilson_prop_field wpcopy, wilson_prop_field wp);
void copy_wv_from_swv(wilson_vector *wv, spin_wilson_vector *swv, int spin);
wilson_prop_field create_wp_field_copy(wilson_prop_field w);
void destroy_swv_field(spin_wilson_vector *swv);
void destroy_wv_field(wilson_vector *wv);
void destroy_wp_field(wilson_prop_field wp);
void destroy_v_field(su3_vector *v);

/* make_prop.c */
int get_wprop_to_wp_field(int startflag, char startfile[], 
			  int saveflag, char savefile[],
			  wilson_prop_field wp,
			  wilson_quark_source *my_wqs,
			  quark_invert_control *my_qic,
			  dirac_clover_param *my_dcp,
			  int check);

int get_ksprop_to_wp_field(int startflag, char startfile[], 
			   int saveflag, char savefile[],
			   wilson_prop_field wp,
			   ks_quark_source *my_ksqs,
			   quark_invert_control *my_qic,
			   ks_param *my_ksp,
			   int check);
/* spectrum_cl.c */
int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
void spectrum_cl(wilson_prop_field qp0, wilson_prop_field qp1, int pair);

/*  cl_inv_includes.h */




