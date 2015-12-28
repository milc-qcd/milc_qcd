/******* ks_multicg_offset_qphix_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QPHIX inverter */

/* C. DeTar 12/26/15 created */

/* NOTE: This code is actually an include file for ks_multicg_qphix_F.c
   and ks_multicg_qphix_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_MULTICG_OFFSET_FIELD

*/

#if ( QPHIX_PrecisionInt == 1 )

#define KS_MULTICG_OFFSET_FIELD ks_multicg_offset_field_qphix_F
#define MYREAL QPHIX_F_Real
#define MYSU3_MATRIX fsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_F_G
#define QPHIX_asqtad_invert QPHIX_F3_asqtad_invert
#define unload_qphix_V_to_field unload_qphix_F_V_to_field

#else

#define KS_MULTICG_OFFSET_FIELD ks_multicg_offset_field_qphix_D
#define MYREAL QPHIX_D_Real
#define MYSU3_MATRIX dsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_D_G
#define QPHIX_asqtad_invert QPHIX_D3_asqtad_invert
#define unload_qphix_V_to_field unload_qphix_D_V_to_field

#endif


#include "generic_ks_includes.h"
#include "../include/generic_qphix.h"
#include "../include/generic_ks_qphix.h"
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <assert.h>

/* Load inversion args for Level 3 inverter */

static void 
set_qphix_invert_arg( QPHIX_invert_arg_t* qphix_invert_arg, 
		    quark_invert_control *qic, int nmass )
{
  qphix_invert_arg->parity       = milc2qphix_parity(qic->parity);
  qphix_invert_arg->max          = qic->max;
  qphix_invert_arg->nrestart     = qic->nrestart;

  /* For multimass inversion, don't restart */
  if(nmass != 1)
    qphix_invert_arg->nrestart = 1;
}


static QPHIX_resid_arg_t **
create_qphix_resid_arg( quark_invert_control *qic, int nmass)
{
  QPHIX_resid_arg_t **res_arg;
  char myname[] = "create_qphix_resid_arg";
  int imass;

  /* Pointers for residual errors */
  res_arg = (QPHIX_resid_arg_t **)malloc(num_offsets*sizeof(QPHIX_resid_arg_t *));
  for(imass = 0; imass < nmass; imass++){
    res_arg[imass] = (QPHIX_resid_arg_t *)malloc(sizeof(QPHIX_resid_arg_t ));
    if(res_arg[imass] == NULL){
      printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
      terminate(1);
    }
    *res_arg[imass] = QPHIX_RESID_ARG_DEFAULT;
  }
  /* For now the residuals are the same for all sources and masses */
  res_arg->resid = qic->resid * qic->resid;
  res_arg->relresid     = 0.;  /* NOT SUPPORTED */
  res_arg->final_rsq    = 0.;
  res_arg->final_rel    = 0.;

  return res_arg;
}

static void
destroy_qphix_resid_arg(QPHIX_resid_arg_t **res_arg, int nmass)
{
  int imass;

  for(imass = 0; imass < nmass[isrc]; imass++){
    free(res_arg[imass]);
  }
  free(res_arg);
}

int
KS_MULTICG_OFFSET_FIELD(
  su3_vector *src,            /* source vector (type su3_vector) */
  su3_vector *psim[],         /* solution vectors */
  ks_param *ksp,	      /* the offsets */
  int nmass,	              /* number of offsets */
  quark_invert_control *qic,/* inversion parameters */
  imp_ferm_links_t *fn      /* Storage for fat and Naik links */
			)
{
  int i,j;
  char myname[] = "ks_multicg_offset_field_qphix";

  // Suggested improvement:
  // Separate the input and output parameters
  // QPHIX_invert_arg_t qphix_invert_arg; (Input)
  // QPHIX_resid_arg_t qphix_resid_arg;  (Output)
  int num_iters = 0; // number of iterations taken
  QPHIX_info_t info = QPHIX_INFO_ZERO;
  MYREAL mass[nmass];
  QPHIX_ColorVector *qphix_sol[nmass], *qphix_src;
  QPHIX_resid_arg_t  **qphix_resid_arg;
  QPHIX_invert_arg_t qphix_invert_arg = QPHIX_INVERT_ARG_DEFAULT;
  QPHIX_FermionLinksAsqtad  *links;    

  /* Initialize QPHIX if not already done */
  if(initialize_qphix() != QPHIX_SUCCESS){
    node0_printf("%s: Error initializing QPHIX\n",myname);
    terminate(1);
  }

#ifdef CGTIME
  double dtimec = -dclock();
  double nflop = 1205 + 15*nmass;
#endif

  assert(qic[0].parity != EVENANDODD && "EVENANDODD not yet implemented");
    
  if(qic[0].relresid != 0.){
    printf("%s: QPhiX code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
  }

  if( nmass==0 )return 0;

  /* Set qphix_invert_arg */
  set_qphix_invert_arg( & qphix_invert_arg, qic+0, nmass );

  /* Pointers for residual errors */
  qphix_resid_arg = create_qphix_resid_arg( nmass, qic+0 );

  /* Map the masses */
  for(i = 0; i < nmass; i++)
    mass[i] = sqrt(ksp[i].offset/4.0);

  /* Map the input and output fields */

  qphix_src = create_V_from_field( src, qic[0].parity);

  
  for(imass = 0; imass < nmass; imass++){
    qphix_sol[imass] = 
      create_V_from_field( psim[imass], qic[0].parity);
  }

  links = create_L_from field( fn, qic->parity );

#ifdef CG_DEBUG
  node0_printf("Calling qphix_ks_multicg_offset\n");fflush(stdout);
#endif

  num_iters = QPHIX_asqtad_multi_invert( info, links, inv_arg, res_arg, mass, 
					 nmass, qphix_dest, qphix_src );

  /* For now we don't support separate residuals for each mass */
  for(i=0; i<nmass; ++i){
    qic[i].final_rsq = resid_arg[i]->final_rsq;
    qic[i].final_relrsq = 0.; /* Not supported at the moment */
    qic[i].final_iters = num_iters;
    qic[i].size_r = resid_arg[i]->size_r;
    qic[i].size_relr = resid_arg[i]->size_relr;
  }

    /* Free the structure */
    destroy_qphix_resid_arg(resid_arg);

  /* Unpack the solutions */
#ifdef CG_DEBUG
  node0_printf("Extracting output\n");fflush(stdout);
#endif

  for(i=0; i<nmass; ++i)
    /* Copy results back to su3_vector */
    unload_qphix_V_to_field( sol[i], qphix_dest[i], qic->parity);
    

    /* Free QPHIX fields  */
    
    QPHIX_destroy_V(qphix_src);    
    QPHIX_destroy_V(qphix_src);    
    for(imass = 0; imass < nmass; imass++)
      QPHIX_destroy_V(qphix_sol[imass]);     

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    char *prec_label[2] = {"F", "D"};
    printf("CONGRAD5: time = %e (multicg_offset_QPHIX %s) masses = %d iters = %d mflops = %e\n",
	   dtimec,prec_label[qic[0].prec-1],num_offsets,num_iters,
	   (double)(nflop)*volume*
	   num_iters/(1.0e6*dtimec*numnodes()));
    fflush(stdout);}
#endif

  return num_iters;
}
