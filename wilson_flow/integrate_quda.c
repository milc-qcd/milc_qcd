/************************* integrate_quda.c **********************/
/* Wilson flow with QUDA */

#include "wilson_flow_includes.h"
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"


void
run_gradient_flow_quda() {

  node0_printf("Running gradient flow with QUDA\n");

  /* Initialize QUDA parameters */
  initialize_quda();

  /* Get gauge field */
  su3_matrix *links = create_G_from_site_quda();

  /* Setup QUDA gauge parameters */
  QudaGaugeParam qgp = newQudaGaugeParam();

  const int * nsquares = get_logical_dimensions();
  qgp.struct_size = sizeof(QudaGaugeParam);
  qgp.type = QUDA_SU3_LINKS;
  qgp.X[0] = nx / nsquares[0];
  qgp.X[1] = ny / nsquares[1];
  qgp.X[2] = nz / nsquares[2];
  qgp.X[3] = nt / nsquares[3];
  qgp.cpu_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qgp.cuda_prec = qgp.cpu_prec;
  qgp.cuda_prec_sloppy = qgp.cuda_prec;
  qgp.cuda_prec_precondition = qgp.cuda_prec;
  qgp.cuda_prec_eigensolver = qgp.cuda_prec;
  qgp.cuda_prec_refinement_sloppy = qgp.cuda_prec;
  qgp.reconstruct = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_eigensolver = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;
  qgp.gauge_order = QUDA_MILC_GAUGE_ORDER;
  qgp.anisotropy = 1.0;
  qgp.t_boundary = QUDA_PERIODIC_T;
  qgp.gauge_fix = QUDA_GAUGE_FIXED_NO;
  qgp.staggered_phase_type = QUDA_STAGGERED_PHASE_NO; // ???i
int pad_size = 0;
  int x_face_size = qgp.X[1] * qgp.X[2] * qgp.X[3] / 2;
  int y_face_size = qgp.X[0] * qgp.X[2] * qgp.X[3] / 2;
  int z_face_size = qgp.X[0] * qgp.X[1] * qgp.X[3] / 2;
  int t_face_size = qgp.X[0] * qgp.X[1] * qgp.X[2] / 2;
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
  pad_size = MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  
  qgp.ga_pad = pad_size; 
  qgp.mom_ga_pad = 0;

  /* Load gauge field in QUDA */
  loadGaugeQuda( (void*) links, &qgp );

  /* Setup QUDA smearing parameters */
  QudaGaugeSmearParam smearParams = newQudaGaugeSmearParam();

  smearParams.struct_size = sizeof(QudaGaugeSmearParam);
  smearParams.n_steps = stoptime / stepsize;
  smearParams.epsilon = stepsize;
  smearParams.meas_interval = 1;
  if( strcmp("wilson", flow_description) == 0 ) {
    smearParams.smear_type = QUDA_GAUGE_SMEAR_WILSON_FLOW;
  }
  else if( strcmp("symanzik", flow_description) == 0 ) {
    smearParams.smear_type = QUDA_GAUGE_SMEAR_SYMANZIK_FLOW;
  }
  else {
    node0_printf("ERROR: flow_description %s is invalid for QUDA flow\n",
           flow_description);
    terminate(1);
  }
  smearParams.restart = QUDA_BOOLEAN_FALSE; // ???
  smearParams.t0 = 0;
  
#ifdef ANISOTROPY
  smearParams.smear_anisotropy = ani;
#endif

#if GF_INTEGRATOR==INTEGRATOR_LUSCHER
  smearParams.rk_order = 3;
#elif GF_INTEGRATOR==INTEGRATOR_BBB
  smearParams.rk_order = 4;
#else
  printf( "run_gradient_flow_quda: unsupported integrator for QUDA\n");
  fflush(stdout); terminate(1);
#endif

  /* Setup QUDA observable parameters */
  int nObsParams = smearParams.n_steps / smearParams.meas_interval + 1;
  QudaGaugeObservableParam *obsParams;
  obsParams = (QudaGaugeObservableParam *)malloc(nObsParams*sizeof(QudaGaugeObservableParam));
  for( int i=0; i<nObsParams; i++)
  {
    obsParams[i] = newQudaGaugeObservableParam();
    obsParams[i].struct_size = sizeof(QudaGaugeObservableParam);
    obsParams[i].su_project = QUDA_BOOLEAN_FALSE;
    obsParams[i].compute_plaquette = QUDA_BOOLEAN_TRUE;
    obsParams[i].compute_rectangle = QUDA_BOOLEAN_TRUE;
    obsParams[i].compute_polyakov_loop = QUDA_BOOLEAN_FALSE;
    obsParams[i].compute_qcharge = QUDA_BOOLEAN_TRUE;
    obsParams[i].compute_qcharge_density = QUDA_BOOLEAN_FALSE;
    obsParams[i].compute_gauge_loop_trace = QUDA_BOOLEAN_FALSE;
  }

  /* Do the gauge flow */
  performWFlowQuda(&smearParams, obsParams);

  /* Clean up */
  destroy_G_quda(links);
  free(obsParams);
}
