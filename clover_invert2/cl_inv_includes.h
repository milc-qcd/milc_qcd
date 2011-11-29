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
#define ENDTIME(string) dtime += dclock(); node0_printf("Aggregate time to %s %e\n",(string),dtime);
#else
#define STARTTIME
#define ENDTIME(string)
#endif

/* prototypes for functions in high level code */
int setup(void);
int readin(int prompt);
int num_mes_report(void);

/* ksprop_info.c */
char *create_ks_XML(void);

/* make_prop.c */
int get_wprop_to_wp_field(int startflag, char startfile[], 
			  int saveflag, char savefile[],
			  wilson_prop_field *wp,
			  quark_source *my_wqs,
			  quark_invert_control *my_qic,
			  dirac_clover_param *my_dcp,
			  Real bdry_phase[],
			  int r0[],
			  int check);

int get_ksprop_to_wp_field(int startflag, char startfile[], 
			   int saveflag, char savefile[],
			   wilson_prop_field *wp,
			   quark_source *my_ksqs,
			   quark_invert_control *my_qic,
			   ks_param *my_ksp,
			   Real *bdry_phase,
			  int r0[],
			   int check);

int get_ksprop4_to_wp_field(int startflag, char startfile[], 
			    int saveflag, char savefile[],
			    wilson_prop_field *wp,
			    quark_source *my_ksqs,
			    quark_invert_control *my_qic,
			    ks_param *my_ksp,
			    Real bdry_phase[],
			    int r0[],
			    int check);

void dump_wprop_from_wp_field(int saveflag, char savefile[], 
			      wilson_prop_field *wp);
void reread_wprop_to_wp_field(int saveflag, char savefile[], wilson_prop_field *wp);

/* spectrum_cl.c */
int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
void spectrum_cl(wilson_prop_field *qp0, wilson_prop_field *qp1, int pair);

/* ks_source_info.c */
char *create_kss_XML(char *filename, quark_source *ksqs);
void free_kss_XML(char *xml);

/* w_source_info.c */
char *create_ws_XML(char *filename, quark_source *wqs);
void free_ws_XML(char *info);

/*  cl_inv_includes.h */




