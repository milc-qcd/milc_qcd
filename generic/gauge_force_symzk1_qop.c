/********************* gauge_force_symzk1_qop.c  -- *************************/
/* MIMD version 7 */
/* Wrapper for QOP gauge force for Symanzik 1 loop gauge action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* Ludmila Levkova created 10/06 
* C.D. Moved to separate file 10/06  */

#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include <qop.h>

void load_qop_imp_gauge_coeffs(QOP_gauge_coeffs_t *c)
{
  Real **loop_coeff = get_loop_coeff();
  int nloop = get_nloop();
  c->plaquette = loop_coeff[0][0];
  if(nloop>1)
    c->rectangle = loop_coeff[1][0];
  else
    /* Hack to support some unimproved actions */
    c->rectangle = 0.;
  if(nloop>2)
    /* Hack to support some unimproved actions */
    c->parallelogram  = loop_coeff[2][0];
  else
    c->parallelogram = 0.;
  //printf("Coefficients loaded %e %e %e\n", c->plaq, c->rect, c->pgm);
}

#ifdef GFTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

void imp_gauge_force( Real eps, field_offset mom_off ){

  QOP_GaugeField *links;
  QOP_Force *mom;

  QOP_gauge_coeffs_t coeff;
  QOP_info_t info = {0., 0., 0, 0, 0};
  double remaptime = -dclock();

  //printf("Begin wrapper\n");
  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("imp_gauge_force: Error initializing QOP\n");
    terminate(1);
  }
  //QDP_Subset *site0_array = QDP_create_subset(S0,NULL,0,1);
  //QDP_Subset site0 = site0_array[0];

 /* Load gauge links and momentum */
  links = create_G_from_site4(F_OFFSET(link), EVENANDODD);
  mom   = create_F_from_site4(F_OFFSET(mom), EVENANDODD);

  /* Load coefficients */
  load_qop_imp_gauge_coeffs(&coeff);
  //	printf("Before calling QOP force\n");

#if 0
  /* Debugging */
  {
#include <qop_internal.h>
    QLA_ColorMatrix *mom0;
    const int x[4] = {0,0,0,0};
    int ind = QDP_index(x);
    
    mom0 = QDP_expose_M(mom->force[0]);        
    node0_printf("FORCE:\n%e %e %e\n\n", mom0[ind].e[0][0].real, mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);
    
    QDP_reset_M(mom->force[0]);
  }
#endif
  remaptime += dclock();
  /* Compute fermion force */
  QOP_symanzik_1loop_gauge_force(&info, links, mom, &coeff, eps*beta);
  remaptime -= dclock();
  //       printf("After calling QOP force status is %d\n",info.status);
  if(info.status != QOP_SUCCESS){
    terminate(1);
  }
  /* Unload momentum */
  unload_F_to_site4( F_OFFSET(mom), mom, EVENANDODD );

  QOP_destroy_G(links); links = NULL;
  QOP_destroy_F(mom);   mom = NULL;

  remaptime += dclock();
#ifdef GFTIME
  node0_printf("GFTIME:  time = %e (qop %s) mflops = %e\n",
	       info.final_sec, qop_prec[QOP_Precision-1],
	       info.final_flop/(1e6*info.final_sec) );
#ifdef REMAP
  node0_printf("GFREMAP:  time = %e\n",remaptime);
#endif
#endif
  //   printf("imp_gauge_force: using my QOP version\n");	
} /* imp_gauge_force.c */


