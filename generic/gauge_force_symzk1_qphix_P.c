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
#include "../include/generic_qphixmilc.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <assert.h>
#include <math.h>

/* Redefinitions according to requested QPhiX precision */

#if ( QPHIX_PrecisionInt == 1 )
#define CREATE_F_FROM_SITE4       create_qphix_F_F_from_site4
#define CREATE_G_FROM_SITE4       create_qphix_F_G_from_site4
#define CREATE_RAW_F_FROM_SITE4   create_qphix_raw4_F_F_from_site4
#define CREATE_RAW_G              create_qphix_raw4_F_G
#define CREATE_RAW_G_FROM_SITE4   create_qphix_raw4_F_G_from_site4
#define DESTROY_RAW_F             destroy_qphix_raw4_F_F
#define DESTROY_RAW_G             destroy_qphix_raw4_F_G
#define IMP_GAUGE_FORCE_QPHIX     imp_gauge_force_qphix_F
#define MYREAL                    float
#define MYSU3_MATRIX              fsu3_matrix
#define UNLOAD_F_TO_SITE4         unload_qphix_F_F_to_site4

#else

#define CREATE_F_FROM_SITE4       create_qphix_D_F_from_site4
#define CREATE_G_FROM_SITE4       create_qphix_D_G_from_site4
#define CREATE_RAW_F_FROM_SITE4   create_qphix_raw4_D_F_from_site4
#define CREATE_RAW_G              create_qphix_raw4_D_G
#define CREATE_RAW_G_FROM_SITE4   create_qphix_raw4_D_G_from_site4
#define DESTROY_RAW_F             destroy_qphix_raw4_D_F
#define DESTROY_RAW_G             destroy_qphix_raw4_D_G
#define IMP_GAUGE_FORCE_QPHIX     imp_gauge_force_qphix_D
#define MYREAL                    double
#define MYSU3_MATRIX              dsu3_matrix
#define UNLOAD_F_TO_SITE4         unload_qphix_D_F_to_site4

#endif

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

/*!
 * Copy backward raw links/momentum without the adjoint.
 * QPhiX does the adjoint in the generated code for the dslash kernels for the back
 * links.
 */
static MYSU3_MATRIX *
create_backlinks_without_adjoint(MYSU3_MATRIX *t)
{
  MYSU3_MATRIX *t_bl = NULL;
  register int i;
  register site *s;
  int dir;
  MYSU3_MATRIX *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "create_backlinks_without_adjoint";
  
  /* Allocate space for t_lbl */
  t_bl = CREATE_RAW_G();
  if(t_bl==NULL){
    printf("%s(%d): no room for t_lbl\n",myname,this_node);
    terminate(1);
  }
  
  tempmat1 = CREATE_RAW_G();
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }
  
  /* gather backwards fatlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLFIELDSITES_OMP(i,){
      tempmat1[i] = t[dir+4*i];
    } END_LOOP_OMP
    tag[dir] = start_gather_field( tempmat1
				   , sizeof(su3_matrix)
				   , OPP_DIR(dir)
				   , EVENANDODD
				   , gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLFIELDSITES_OMP(i,) {
      MYSU3_MATRIX * temp_ = (t_bl + dir + 4*i);
      t_bl[dir + 4*i] = *((MYSU3_MATRIX *)gen_pt[dir][i]);
    } END_LOOP_OMP
    cleanup_gather( tag[dir] );
  }
  
  DESTROY_RAW_G(tempmat1); 
  tempmat1 = NULL;
  return t_bl;
}

void
IMP_GAUGE_FORCE_QPHIX ( Real eps, field_offset mom_off )
{
  char myname[] = "imp_gauge_force_qphix";
  int nflop = 153004;  /* For Symanzik1 action, based on MILC standard */
  QPHIX_Force *mom;
  MYSU3_MATRIX *fwdrawmom;
  QPHIX_GaugeField *gauge;
  MYSU3_MATRIX *fwdrawgauge;
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
  fwdrawmom = CREATE_RAW_F_FROM_SITE4(F_OFFSET(mom), EVENANDODD);
  //  bckrawmom = create_backlinks_without_adjoint(fwdrawmom);
  //  mom = QPHIX_create_F_from_raw((MYREAL *)fwdrawmom, (MYREAL *)bckrawmom, milc2qphix_parity(EVENANDODD));
  mom = QPHIX_create_F_from_raw((MYREAL *)fwdrawmom, milc2qphix_parity(EVENANDODD));
  DESTROY_RAW_F(fwdrawmom);
  //  DESTROY_RAW_F(bckrawmom);

  /* Convert the gauge field to QPhiX format */
  fwdrawgauge = CREATE_RAW_G_FROM_SITE4(F_OFFSET(link[0]), EVENANDODD);
  //  bckrawgauge = create_backlinks_without_adjoint(fwdrawgauge);
  //  gauge = QPHIX_create_G_from_raw((MYREAL *)fwdrawgauge, (MYREAL *)bckrawgauge, milc2qphix_parity(EVENANDODD));
  gauge = QPHIX_create_G_from_raw((MYREAL *)fwdrawgauge, milc2qphix_parity(EVENANDODD));
  DESTROY_RAW_G(fwdrawgauge);
  //  DESTROY_RAW_G(bckrawgauge);

  /* Update the mom, based on the gauge force */
  QPHIX_symanzik_1loop_gauge_force( &info, gauge, mom, coeffs, eps*beta/3.);

  if(info.status != 0){
    node0_printf("Quitting because of gauge force error\n");
    terminate(1);
  }

  /* Copy the result to the site structure */

  UNLOAD_F_TO_SITE4(F_OFFSET(mom), mom, EVENANDODD);
  
  /* Clean up */

  QPHIX_destroy_F(mom);
  QPHIX_destroy_G(gauge);

#ifdef GFTIME
  dtime+=dclock();
  qtime = info.final_sec;
  node0_printf("GFTIME_QPHIX:   time = %e (QPhiX Symanzik1) mflops = %e\n",
	       qtime, 
	       info.final_flop*(double)volume/(1e6*qtime*numnodes()) );
  node0_printf("GFTIME:   time = %e (QPhiX Symanzik1) mflops = %e\n",
	       dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif
}
