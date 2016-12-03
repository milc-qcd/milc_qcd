/******* ks_multicg_offset_grid_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 GRID inverter */

/* C. DeTar 12/26/15 created */

/* NOTE: This code is actually an include file for ks_multicg_grid_F.c
   and ks_multicg_grid_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_MULTICG_OFFSET_FIELD

*/

#if ( GRID_PrecisionInt == 1 )

#define KS_MULTICG_OFFSET_FIELD ks_multicg_offset_field_grid_F
#define MYREAL GRID_F_Real
#define MYSU3_MATRIX fsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_F_G
#define unload_grid_V_to_field unload_grid_F_V_to_field

#else

#define KS_MULTICG_OFFSET_FIELD ks_multicg_offset_field_grid_D
#define MYREAL GRID_D_Real
#define MYSU3_MATRIX dsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_D_G
#define unload_grid_V_to_field unload_grid_D_V_to_field

#endif


#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <assert.h>
#ifdef VTUNE
#include <ittnotify.h>
#endif

/* Load inversion args for Level 3 inverter */

static void 
set_grid_invert_arg( GRID_invert_arg_t* grid_invert_arg, 
		    quark_invert_control *qic, int nmass )
{
  grid_invert_arg->parity       = milc2grid_parity(qic->parity);
  grid_invert_arg->max          = qic->max;
  grid_invert_arg->nrestart     = qic->nrestart;

  /* For multimass inversion, don't restart */
  if(nmass != 1)
    grid_invert_arg->nrestart = 1;
}


static GRID_resid_arg_t **
create_grid_resid_arg( quark_invert_control *qic, int nmass)
{
  GRID_resid_arg_t **res_arg;
  char myname[] = "create_grid_resid_arg";
  int imass;

  /* Pointers for residual errors */
  res_arg = (GRID_resid_arg_t **)malloc(nmass*sizeof(GRID_resid_arg_t *));
  for(imass = 0; imass < nmass; imass++){
    res_arg[imass] = (GRID_resid_arg_t *)malloc(sizeof(GRID_resid_arg_t ));
    if(res_arg[imass] == NULL){
      printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
      terminate(1);
    }
    *res_arg[imass] = GRID_RESID_ARG_DEFAULT;
  }
  /* For now the residuals are the same for all sources and masses */
  for(imass = 0; imass < nmass; imass++){
    res_arg[imass]->resid = qic->resid * qic->resid;
    res_arg[imass]->relresid     = 0.;  /* NOT SUPPORTED */
    res_arg[imass]->final_rsq    = 0.;
    res_arg[imass]->final_rel    = 0.;
  }

  return res_arg;
}

/* Collect inversion statistics */

static void 
get_grid_resid_arg( quark_invert_control *qic, 
		     GRID_resid_arg_t **grid_resid_arg, int nmass, int num_iters )
{
  /* For now we don't support separate residuals for each mass */
  for(int i=0; i<nmass; ++i){
    qic[i].final_rsq     = grid_resid_arg[i]->final_rsq;
    qic[i].final_relrsq  = 0.;                            /* Not supported at the moment */
    qic[i].final_iters   = num_iters;
    qic[i].size_r        = grid_resid_arg[i]->size_r;
    qic[i].size_relr     = grid_resid_arg[i]->size_relr;
    qic[i].final_iters   = grid_resid_arg[i]->final_iter;
    qic[i].final_restart = grid_resid_arg[i]->final_restart;
  }
}

static void
destroy_grid_resid_arg(GRID_resid_arg_t **res_arg, int nmass)
{
  int imass;

  for(imass = 0; imass < nmass; imass++){
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
  char myname[] = "ks_multicg_offset_field_grid";

  // Suggested improvement:
  // Separate the input and output parameters
  // GRID_invert_arg_t grid_invert_arg; (Input)
  // GRID_resid_arg_t grid_resid_arg;  (Output)
  int num_iters = 0; // number of iterations taken
  GRID_info_t info = GRID_INFO_ZERO;
  MYREAL mass[nmass];
  GRID_ColorVector *grid_sol[nmass], *grid_src;
  GRID_resid_arg_t  **grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_FermionLinksAsqtad  *links;    

  /* Initialize GRID if not already done */
  if(initialize_grid(GRID_PrecisionInt) != GRID_SUCCESS){
    node0_printf("%s: Error initializing GRID\n",myname);
    terminate(1);
  }

#ifdef CGTIME
  double dtimec = -dclock();
  double nflop = 1205 + 15*nmass;
#endif

  assert(qic[0].parity != EVENANDODD && "EVENANDODD not yet implemented");
    
  if(qic[0].relresid != 0.){
    node0_printf("%s: Grid code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
  }

  if( nmass==0 )return 0;

#ifdef VTUNE
  node0_printf("Starting VTune\n");
  __itt_resume();
#endif
  double dtimet = -dclock();
  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic+0, nmass );
  node0_printf("set_grid_invert_arg: time = %.6e sec\n", dclock()+dtimet);

  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic+0, nmass );

  /* Map the masses */
  for(i = 0; i < nmass; i++)
    mass[i] = sqrt(ksp[i].offset/4.0);

  /* Map the input and output fields */
  dtimet = -dclock();
  grid_src = create_grid_V_from_field( src, qic[0].parity);
  node0_printf("create_grid_V_from_field: time = %.6e sec\n", dclock()+dtimet);
  fflush(stdout);
  
  dtimet = -dclock();
  for(i = 0; i < nmass; i++){
    grid_sol[i] = 
      create_grid_V_from_field( psim[i], qic[0].parity);
  }
  node0_printf("create_grid_V_from_field x %d: time = %.6e sec\n", 
	       nmass, dclock()+dtimet);
  fflush(stdout);

  dtimet = -dclock();
  links = create_grid_L_from_fn_links( fn, EVENANDODD );
  node0_printf("create_grid_L_from_fn_links: time = %.6e sec\n", dclock()+dtimet);
  fflush(stdout);

#ifdef CG_DEBUG
  node0_printf("Calling GRID_ks_multicg_offset\n");fflush(stdout);
#endif

  num_iters = GRID_asqtad_invert_multi( &info, links, &grid_invert_arg, grid_resid_arg, 
					 mass, nmass, grid_sol, grid_src );

  get_grid_resid_arg( qic, grid_resid_arg, nmass, num_iters);

  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg, nmass);
  
  /* Unpack the solutions */
#ifdef CG_DEBUG
  node0_printf("Extracting output\n");fflush(stdout);
#endif

  dtimet = -dclock();
  for(i=0; i<nmass; ++i)
    /* Copy results back to su3_vector */
    unload_grid_V_to_field( psim[i], grid_sol[i], qic->parity);
  node0_printf("unload_grid_V_to_field x %d: time = %.6e sec\n", 
	       nmass, dclock()+dtimet);
  fflush(stdout);
  
  /* Free GRID fields  */
  
  GRID_destroy_V(grid_src);    
  for(i = 0; i < nmass; i++)
    GRID_destroy_V(grid_sol[i]);     

  GRID_asqtad_destroy_L(links);
  
#ifdef VTUNE
  __itt_pause();
  node0_printf("Ending VTune\n");
#endif

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    char *prec_label[2] = {"F", "D"};
    printf("CONGRAD5: time = %e (multicg_offset_GRID %s) masses = %d iters = %d mflops = %e\n",
	   dtimec,prec_label[qic[0].prec-1],nmass,num_iters,
	   (double)(nflop)*volume*
	   num_iters/(1.0e6*dtimec*numnodes()));
    fflush(stdout);}
#endif

  total_iters += num_iters;
  return num_iters;
}
