/****************** wi_dyn_includes.h ******************************/
/*
*  Include files for the wilson_invert application
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
#include "../include/io_lat.h"
#include "../include/io_wb.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"

/* prototypes for functions in high level code */
int setup_w();
int readin(int prompt);
/*  wi_inv_includes.h */

