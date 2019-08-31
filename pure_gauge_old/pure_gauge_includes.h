/****************** pure_gauge_includes.h ******************************/
/*
*  Include files for pure gauge (Wilson action) application
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
#include "../include/generic_pg.h"
#include "../include/dirs.h"
#include <string.h>

/* prototypes for functions in high level code */
int  setup();
int readin(int prompt);

