/****************** string_break_includes.h ******************************/
/*
*  Include files for heavy quark potential application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/dirs.h"

int setup();
void smearing( void );
int readin(int prompt);
void ax_gauge();
void fuz_source(field_offset grd, field_offset fsrc, int r0, int t0, int step);
void fuz_prop(field_offset fprop, int r0);
void stat_li_mesons(int tot_smear, int step);
void li_li_mesons(int tot_smear, int t0, char *sink);
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
void wl_1l_1corr(int tot_smear, int step);
void wl_1l_2corr(int tot_smear, int step);
void wl_2l_1corr(int tot_smear, int step);
void wl_2l_2corr(int tot_smear, int step);
void wl_1l_1corr_offax(int tot_smear, int step);
void wl_1l_2corr_offax(int tot_smear, int step);
void wl_2l_1corr_offax(int tot_smear, int step);
void wl_2l_2corr_offax(int tot_smear, int step);

void su3vecsrc_copy(su3_vector_src *a, su3_vector_src *b, int num_src);
void su3_vec_to_src(su3_vector *a, su3_vector_src *b, int num_src);
void su3_src_to_vec(su3_vector_src *a, su3_vector *b, int j);
void su3vecsrc_outer_prod(su3_vector *a, su3_vector_src *b,
			  su3_matrix *c, int num_src);

