/****************** schroed_ks_includes.h ******************************/
/*
*  Include files for Kogut-Susskind Schroedinger functional
*  dynamical application
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
#include "../include/generic_ks.h"
#include "../include/generic_schroed.h"
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int  setup();
int readin(int prompt);
int update();
void update_h(Real eps);
void update_u(Real eps);
void gauge_force(Real eps);
void fermion_force(Real eps);
double d_action();

void gauge_field_copy(field_offset src, field_offset dest);

void f_measure(Real *r_psi_bar_psi_even, Real *r_psi_bar_psi_odd,
	       Real *r_ferm_energy, Real *r_ferm_pressure,
	       Real *r_ferm_action);

void phaseset_sf();
void rephase_sf(int flag);
