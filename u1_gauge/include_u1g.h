/* ************************************************************ */
/*								*/
/*				INCLUDE_U1G.H			*/
/*								*/
/* To be included in most of the generic routine files		*/
/*								*/
/* Last Updated on 04.18.07					*/
/*								*/
/* ************************************************************ */

#include "../include/config.h"  /* keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_u1lat.h"
#include "../include/generic.h"
#include "../include/generic_u1.h"
#include "../include/random.h"
#include "../include/dirs.h"

#include "lattice.h"

/* Compile-time defaults for this application.  The field generator in
   generic_u1/u1_gauge_field.c takes these as run-time arguments, so an
   external caller need not define them. */
#ifndef GAUGE
#define GAUGE U1_COULOMB
#endif
#ifndef SCHEME
#define SCHEME U1_QED_L
#endif

int setup(void);
int readin(int prompt);

/* ************************************************************ */

