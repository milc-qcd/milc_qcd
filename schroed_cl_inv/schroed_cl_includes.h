/****************** schroed_cl_includes.h ******************************/
/*
*  Include files for clover Schroedinger functional inverter application
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
#include "../include/generic_schroed.h"
#include "../include/dirs.h"
#include "../include/check_malloc.h"

/* prototypes for functions in high level code */
int setup_cl();
int readin(int prompt);

void do_phases();
void schroed_meson(field_offset src1, field_offset src2,
		   complex *prop[10], int max_prop);

void ape_smear_SF();
void zv_meas(field_offset src1, field_offset src2,
	     Real *f_1, complex *f_V, Real kappa);

void w_source_sf_site(field_offset src, my_quark_source *wqs);

