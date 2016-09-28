/******************* gauge_force_imp_qphix_P.c ************************/
/* For the QPhiX interface */
/* MIMD version 7 */

/* 9/27/16 C. DeTar

/* NOTE: This code is actually an include file for gauge_force_imp_qphix_F.c
   and gauge_force_imp_qphix_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   GAUGE_FORCE_IMP_QPHIX

*/

#include "../include/generic_qphix.h"
#include "../include/generic_ks_qphix.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* Redefinitions according to requested QPhiX precision */

#if ( QPHIX_PrecisionInt == 1 )
#define CREATE_F_FROM_SITE4       create_qphix_F_F_from_site4
#define CREATE_G_FROM_SITE4       create_qphix_F_G_from_site4
#define IMP_GAUGE_FORCE_QPHIX     imp_gauge_force_qphix_F
#define UNLOAD_F_TO_SITE4         unload_qphix_F_F_to_site4

#else

#define CREATE_F_FROM_SITE4       create_qphix_D_F_from_site4
#define CREATE_G_FROM_SITE4       create_qphix_D_G_from_site4
#define IMP_GAUGE_FORCE_QPHIX     imp_gauge_force_qphix_D
#define UNLOAD_F_TO_SITE4         unload_qphix_D_F_to_site4

#endif

/* Redefinition according to the prevailing precision */

#if (PRECISION==1)

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

#else

/* Convert (or copy) su3_matrix from prevailing to single precision */
static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

#endif

#define copy_milc_to_F_G(d,s) p2f_mat(d,s);
#define copy_milc_to_D_G(d,s) p2d_mat(d,s);


/* Standard hack to extract macros defined in the action header, */

#define GAUGE_ACTION_PART1
/* defines all loops and their coefficients */
#include <gauge_action.h>
/* The action must be of the symanzik_1loop action type */
#ifndef SYMANZIK_1LOOP
#error "QPhiX supports only Symanzik 1-loop gauge actions"
#endif

#undef GAUGE_ACTION_PART1

/* Load QPHIX_gauge_coeffs from the action header */

static QPHIX_gauge_coeffs_t *
make_gauge_action_coeffs(void){
  static QPHIX_gauge_coeffs_t coeffs;

  /* Symanzik 1-loop actions have only three loops */
  int i,j;
  int total_dyn_flavors;
  char gauge_action_description[128];
  double loop_coeff[3][1];
  int loop_num[3], loop_length[3];
  assert(NLOOP == 3 && NREPS == 1);

  total_dyn_flavors = 0;
  for(i = 0; i < n_dyn_masses; i++){
    total_dyn_flavors += dyn_flavors[i];
  }

  /* The standard ugly hack to get the loop coefficients from the action header */
  /* Defines loop_ind and loop_length_in */
#define GAUGE_ACTION_PART2
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART2

  /* Populate the coeffs structure */
  
  coeffs.plaquette = loop_coeff[0][0];
  coeffs.rectangle = loop_coeff[1][0];
  coeffs.parallelogram = loop_coeff[2][0]; 
  coeffs.adjoint_plaquette = 0.;  /* Not used for Symanzik 1-loop */
    
  return &coeffs;
}

void
IMP_GAUGE_FORCE_QPHIX ( Real eps, field_offset mom_off )
{
  char myname[] = "imp_gauge_force_qphix";
  int nflop = 153004;  /* For Symanzik1 action, based on MILC standard */
  QPHIX_Force *mom;
  QPHIX_GaugeField *gauge;
  QPHIX_info_t info;
  QPHIX_gauge_coeffs_t* coeffs;
#ifdef GFTIME
  double dtime, qtime;
#endif
#ifdef GFTIME
  dtime=-dclock();
#endif

  /* Get gauge action coefficients */
  coeffs = make_gauge_action_coeffs();

  /* Convert the site momentum to QPhiX format */
  mom = CREATE_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);

  /* Convert the gauge field to QPhiX format */
  gauge = CREATE_G_FROM_SITE4(F_OFFSET(link[0]), EVENANDODD);

  /* Update the mom, based on the gauge force */
  QPHIX_symanzik_1loop_gauge_force( &info, gauge, mom, coeffs, eps);

  if(info.status != 0){
    node0_printf("Quitting because of gauge force error\n");
    terminate(1);
  }

  /* Copy the result to the site structure */

  UNLOAD_F_TO_SITE4(F_OFFSET(mom), mom, EVENANDODD);
  
  /* Clean up */

  QPHIX_destroy_F(mom);

#ifdef GFTIME
  dtime+=dclock();
  qtime = info->final_sec;
  node0_printf("GFTIME_QPHIX:   time = %e (QPhiX Symanzik1) mflops = %e\n",
	       qtime, 
	       info->final_flop*(double)volume/(1e6*qtime*numnodes()) );
  node0_printf("GFTIME:   time = %e (QPhiX Symanzik1) mflops = %e\n",
	       dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif
}
