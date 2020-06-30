/*************** gauge_utilities_includes.h ******************************/
/* Include files and prototypes common to most files in the application  */

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

/* Prototypes for functions in high level code */
int setup();
int readin(int prompt);
void flow_step();
//void stout_step_rk();
void staple();
void fmunu_fmunu(double *time, double *space, double *charge);
void initialize_integrator();
//void integrate_RK_2N();
//void integrate_RK_2N_one_step( Real cA, Real cB );
void gauge_action_w_s( double *gact_w_s, double *gact_w_t,
                       double *gact_s_s, double *gact_s_t );
// various integrators, compile-time choice
void integrate_RK_2N();
void integrate_RKMK3();
void integrate_RKMK_generic();
void integrate_adapt_RK_2N();
void integrate_adapt_bs();
