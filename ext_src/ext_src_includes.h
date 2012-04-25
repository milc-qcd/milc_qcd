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
//#include "../include/generic_ks.h"
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
void extract_wprop_to_w_source(int startflag, char startfile[], int ncolor,
			       int nt0, quark_source *dst_qs,
			       quark_source_sink_op *snk_qs_op, int snk_gam,
			       int dst_type);

void extract_ksprop_to_ks_source(int startflag, char startfile[], int ncolor,
				 int nt0, quark_source *dst_qs,
				 quark_source_sink_op *snk_qs_op, 
				 int snk_spin_taste, int r0[]);
void extract_ksprop_to_w_source(int startflag, char startfile[], int ncolor,
				int nt0, quark_source *dst_qs,
				quark_source_sink_op *snk_qs_op, int snk_gam,
				int dst_type);

/* ks_source_info.c */
char *create_kss_XML(char *filename, quark_source *ksqs);
void free_kss_XML(char *xml);

/* w_source_info.c */
char *create_ws_XML(char *filename, quark_source *wqs);
void free_ws_XML(char *info);

/*  ext_src_includes.h */




