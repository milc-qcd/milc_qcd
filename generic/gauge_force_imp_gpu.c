/*********************** gauge_force_imp.c  -- ****************************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* C.D. split from gauge_stuff.c 10/06 */

/**#define GFTIME**/ /* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/generic_quda.h"

// gpu code 
void imp_gauge_force_gpu(Real eps, field_offset mom_off)
{

#ifdef GFTIME
  int nflop = 153004;  /* For Symanzik1 action */
  double dtime = -dclock();
#endif

  Real **loop_coeff = get_loop_coeff();
  //int max_length = get_max_length();
  //int nreps = get_nreps();

  const int num_loop_types = get_nloop();
  double *quda_loop_coeff = (double*)malloc(num_loop_types * sizeof(double));
  int i;
  site *st;
  const Real eb3 = eps*beta/3.0;
  
  initialize_quda();

  su3_matrix *links = create_G_from_site_quda();
  anti_hermitmat* momentum = create_M_quda();

  for(i=0; i<num_loop_types; ++i) quda_loop_coeff[i] = loop_coeff[i][0];

  qudaGaugeForce(PRECISION,num_loop_types,quda_loop_coeff,eb3,links,momentum);

  FORALLSITES_OMP(i,st,){
    for(int dir=XUP; dir<=TUP; ++dir){
      for(int j=0; j<10; ++j){
	((Real*)&(st->mom[dir]))[j] += ((Real*)(momentum + 4*i+dir))[j];
      }
    }
  } END_LOOP_OMP
  
  free(quda_loop_coeff);
  destroy_G_quda(links);
  destroy_M_quda(momentum);

#ifdef GFTIME
  dtime+=dclock();
  node0_printf("GFTIME:   time = %e (Symanzik1_QUDA) mflops = %e\n",dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif

  return;
}

