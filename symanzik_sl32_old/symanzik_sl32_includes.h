/****************** symanzik_sl32_includes.h ******************************/
/*
*  Include files for pure gauge improved actions (needing 32 sublattices)
*  application
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
#include <string.h>
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int  setup();
int readin(int prompt);
int update();
void update_h(Real eps);
void update_u(Real eps);
void relax(int NumStp);
void monte(int NumStp);
double d_action();
double imp_gauge_action();
double hmom_action();
void make_loop_table();
void gauge_field_copy(field_offset src, field_offset dest);

