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

  for(i=0; i<num_loop_types; ++i) quda_loop_coeff[i] = loop_coeff[i][0];

  QudaMILCSiteArg_t arg = newQudaMILCSiteArg();
  qudaGaugeForcePhased(MILC_PRECISION,num_loop_types,quda_loop_coeff,eb3,&arg, phases_in);

  free(quda_loop_coeff);

#ifdef GFTIME
  dtime+=dclock();
  node0_printf("GFTIME:   time = %e (Symanzik1_QUDA) mflops = %e\n",dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif

  return;
}

