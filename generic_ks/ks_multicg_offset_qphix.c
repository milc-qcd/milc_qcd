/************** ks_multicg_offset_qphix.c **************************/
/* MIMD version 7 */

/* The following headers are supplied with the MILC code */
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qphix.h"
#include "../include/loopend.h"
#include <assert.h>

/* qphix multi-mass cg inverter */
/* mbench/ks_multicg_offset.cpp */

// Suggested improved API
//void QPHIX_F3_asqtad_invert(QPHIX_info_t *info,
//			  QPHIX_F3_FermionLinksAsqtad *asqtad,
//			  QPHIX_invert_arg_t *inv_arg,
//			  QPHIX_resid_arg_t *res_arg,
//			  QPHIX_F_Real mass,
//			  QPHIX_F3_ColorVector *out_pt,
//			  QPHIX_F3_ColorVector *in_pt);
//
int qphix_ks_multicg_offset(
			    void *t_src_arg,
			    void *t_dest_arg[],
			    struct QuarkInvertControl *qic,     // as defined in ks_d_congrad_fn.h
			    fptype *mass,
			    int num_offsets,
			    void *gll_arg[2],
			    void *gfl_arg[2]
			    );

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

int ks_multicg_offset_field_qphix(
    su3_vector *src,
    su3_vector **psim,
    ks_param *ksp,
    int num_offsets,
    quark_invert_control *qic,
    imp_ferm_links_t *fn
    )
{
  int i,j;
  char myname[] = "ks_multicg_offset_field_qphix";

  // Suggested improvement:
  // Separate the input and output parameters
  // QPHIX_invert_arg_t qphix_invert_arg; (Input)
  // QPHIX_resid_arg_t qphix_resid_arg;  (Output)
  struct QuarkInvertControl qphix_qic;
  int num_iters = 0; // number of iterations taken

#ifdef CGTIME
  double dtimec = -dclock();
  double nflop = 1205 + 15*num_offsets;
#endif
  int parity = qic[0].parity;

  assert(parity != EVENANDODD && "EVENANDODD not yet implemented");
    
  if(qic[0].relresid != 0.){
    printf("%s: QPhiX code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
  }

  /* Initialize structure */
  for(j = 0; j < num_offsets; j++){
    qic[j].final_rsq     = 0.;
    qic[j].final_relrsq  = 0.; /* No relative residual in use here */
    qic[j].size_r        = 0.;
    qic[j].size_relr     = 0.;
    qic[j].final_iters   = 0;
    qic[j].final_restart = 0;  /* No restarts with this algorithm */
    qic[j].converged     = 1;
  }

  if( num_offsets==0 )return 0;

  /* Remap inversion control structure */
  /* At the moment we are not supporting different residuals for different shifts */
  // Split out into a separate procedure: 
  // set_qphix_invert_arg( QPHIX_invert_arg_t*, qphix_invert_arg, quark_invert_control *qic, int nsrc, int nmass[] );
  qphix_qic.prec      = qic[0].prec;       /* Currently ignored */
  qphix_qic.parity    = parity;
  qphix_qic.max       = qic[0].max;
  qphix_qic.nrestart  = qic[0].nrestart;
  qphix_qic.resid     = qic[0].resid * qic[0].resid;
  qphix_qic.relresid  = qic[0].relresid * qic[0].relresid;          /* Suppresses this test */
  
  /* Map the masses */
  // Replace fptype with QPHIX_F_Real
  fptype* mass = (double*)malloc(num_offsets*sizeof(fptype));
  for(i = 0; i < num_offsets; i++)
    mass[i] = sqrt(ksp[i].offset/4.0);

  /* Map the input and output fields */
  // This allocations should be replaced by a single constructor
  // Suggested constructor 
  // QPHIX_F3_ColorVector *QPHIX_F3_create_V_from_raw( QPHIX_F_Real *src, QPHIX_evenodd_t evenodd);
  // So...
  // QPHIX_F3_ColorVector *qphix_src = QPHIX_F3_create_V_from_raw( src_raw, qphix_parity); 
  void *t_src_arg = (void *)allocKS();
  get_ks_spinors_from_lattice(src, t_src_arg, parity);

  // These allocations should be replaced by a single constructor
  // Suggested constructor
  // QPHIX_F3_FermionLinksAsqtad *
  //  QPHIX_F3_asqtad_create_L_from_raw( QPHIX_F_Real *fatlinks[],
  //			                 QPHIX_F_Real *longlinks[],
  //			                 QPHIX_evenodd_t evenodd);
  // QPHIX_F3_FermionLinksAsqtad *qphix_links = 
  //    QPHIX_F3_asqtad_create_L_from_raw( fatlinks_raw, longlinks_raw, qphix_parity);

  void *gfl_arg[2], *gll_arg[2];
  gfl_arg[0] = allocGauge18();
  get_fatlinks_from_lattice(gfl_arg[0], EVEN, fn);
  gfl_arg[1] = allocGauge18();
  get_fatlinks_from_lattice(gfl_arg[1], ODD, fn);
  gll_arg[0] = allocGauge();
  get_longlinks_from_lattice(gll_arg[0], EVEN, fn);
  gll_arg[1] = allocGauge();
  get_longlinks_from_lattice(gll_arg[1], ODD, fn);
  
  void **t_dest_arg = (void **)malloc(num_offsets*sizeof(void *));
  for(i=0; i<num_offsets; ++i)
    // This allocations should be replaced by a single constructor
    // QPHIX_F3_ColorVector *qphix_sol[i] = QPHIX_F3_create_V_from_raw( sol_raw[i], qphix_parity); 
    t_dest_arg[i] = (void *)allocKS();

  /* Check if the mbench object has been created */
  initialize_qphix();

  num_iters = qphix_ks_multicg_offset(
				      t_src_arg,
				      t_dest_arg,
				      &qphix_qic,     // as defined in ks_d_congrad_fn.h
				      mass,
				      num_offsets,
				      gll_arg,
				      gfl_arg);

  /* Unpack the solutions */
  for(i=0; i<num_offsets; ++i)
    // Suggested accessor
    // void QPHIX_F3_extract_V_to_raw(QPHIX_F_Real *dest, QPHIX_F3_ColorVector *src, QPHIX_evenodd_t evenodd);
    set_ks_spinors_into_lattice (psim[i], t_dest_arg[i], parity);

  /* Unpack the inverter statistics */
  /* For now we don't support separate residuals for each mass */
  for(i=0; i<num_offsets; ++i){
    qic[i].final_rsq = qphix_qic.final_rsq;
    qic[i].final_relrsq = 0.; /* Not supported at the moment */
    qic[i].final_iters = num_iters;
    qic[i].size_r = qphix_qic.size_r;
    qic[i].size_relr = qphix_qic.size_relr;
  }

  free(mass);

  // Suggested destructors
  // void QOP_F3_destroy_V(QOP_F3_ColorVector *field);
  // void QOP_F3_asqtad_destroy_L(QOP_F3_FermionLinksAsqtad *field);

  for(i=0; i<num_offsets; ++i)
    freeKS(t_dest_arg[i]);
  free(t_dest_arg);

  freeGauge(gll_arg[0]);
  freeGauge(gll_arg[1]);
  freeGauge18(gfl_arg[0]);
  freeGauge18(gfl_arg[1]);

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    printf("CONGRAD5: time = %e (multicg_offset_QUDA %s) masses = %d iters = %d mflops = %e\n",
	   dtimec,prec_label[qic[0].prec-1],num_offsets,num_iters,
	   (double)(nflop)*volume*
	   num_iters/(1.0e6*dtimec*numnodes()));
    fflush(stdout);}
#endif

  return num_iters;
}
