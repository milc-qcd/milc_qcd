/****************** ext_src_includes.h ******************************/
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
#include "../include/generic_ks.h"
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

/* make_ext_src.c */
void extract_wprop_to_w_source(int startflag, char startfile[], 
			       int nt0, wilson_quark_source *my_wqs,
			       wilson_quark_source *snk_wqs,
			       int snk_gam);
void extract_ksprop_to_ks_source(int startflag, char startfile[], 
				 int nt0, ks_quark_source *my_ksqs,
				 ks_quark_source *snk_ksqs);
void extract_ksprop_to_w_source(int startflag, char startfile[], 
				int nt0, wilson_quark_source *my_wqs,
				wilson_quark_source *snk_ksqs,
				int snk_gam);
/* ext_src.c */
int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
void spectrum_cl(wilson_prop_field qp0, wilson_prop_field qp1, int pair);

/* ks_source_info.c */
char *create_kss_XML(char *filename, ks_quark_source *ksqs);
void free_kss_XML(char *xml);

/* w_source_info.c */
char *create_ws_XML(char *filename, wilson_quark_source *wqs);
void free_ws_XML(char *info);

/*  ext_src_includes.h */




