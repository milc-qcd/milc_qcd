/****************** hvy_qpot_includes.h ******************************/
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
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

void ax_gauge();
int setup();
int readin(int prompt);
void gball_simp(int tot_smear);
void smearing(void);
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
void hybrid_loop1(int tot_smear);

