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

/* prototypes for functions in high level code */

int  setup();
int readin(int prompt);
void rotate_w_quark(field_offset src, field_offset dest, Real d1);
//void All_KS_hl_prop(field_offset snk, field_offset src, complex **propagator);
FILE* open_fnal_meson_file(char filename[]);
void close_fnal_meson_file(FILE *fp);
void spectrum_hl_rot(FILE *fp, field_offset snk, field_offset src, int k);
void spectrum_hl_smear(FILE *fp, field_offset snk, field_offset src, 
		       int k, int ns);
void get_smearings_bi_serial(char *filename);

/* For baryon_twopt.c */
int ks_baryon_2point_Omu_HH(field_offset ks_1, field_offset heavy_quark, double_complex *propagator[4][4][3][3]);
int ks_baryon_2point_Omu(field_offset ks_1, field_offset ks_2, field_offset heavy_quark, double_complex *propagator[4][4][3][3]);
int ks_baryon_2point_O5(field_offset ks_1, field_offset ks_2, field_offset heavy_quark, double_complex *propagator[4][4]);

