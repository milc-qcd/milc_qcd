/****************** cl_inv_includes.h ******************************/
/*
*  Include files for the clover_invert application
*/

#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
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
void dump_wprop_from_wp_field(int saveflag, char savefile[], 
			      wilson_prop_field wp);
wilson_prop_field reread_wprop_to_wp_field(int saveflag, char savefile[]);

/* spectrum_cl.c */
int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
void spectrum_cl(wilson_prop_field qp0, wilson_prop_field qp1, int pair);

/*  cl_inv_includes.h */




