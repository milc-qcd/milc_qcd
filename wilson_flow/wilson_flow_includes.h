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
void run_gradient_flow();
void flow_step();

//void stout_step_rk();
void staple();
void fmunu_fmunu(double *time, double *space, double *charge);
void initialize_integrator();
void gauge_action_w_s( double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );
// various integrators, compile-time choice
void integrate_RK_2N();
void integrate_RKMK3();
void integrate_RKMK_generic();
void integrate_adapt_RK_2N();
void integrate_adapt_bs();
