/****************** ks_hl_spectrum_includes.h ******************************/
/*
 *  Include files for heavy clover - staggered light spectrum code
 */

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/io_ksprop.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/generic_wilson.h"
#include "../include/dirs.h"
#ifdef HAVE_QDP
#include <qdp.h>
#endif

/* prototypes for functions in high level code */

void weyl2canopy_w_rot(field_offset src, field_offset dest);
void bj_to_w_rot(field_offset src, field_offset dest);
void canopy2weyl_w_rot(field_offset src, field_offset dest);
int  setup();
int readin(int prompt);
void rotate_w_quark(field_offset src, field_offset dest, double d1);
complex  KS_2pt_trace(su3_matrix * antiquark, wilson_propagator * quark, 
		      int * g_snk, int n_snk, int *g_src, int n_src, int *p, site *s);
void KS_2pt_func(field_offset snk, field_offset src, int *g_snk, int n_snk,
		 int *g_src, int n_src, int *p, complex *prop, int parity);
void All_KS_hl_prop(field_offset snk, field_offset src, complex **propagator);
void get_smearings_bi_serial(char *filename);
