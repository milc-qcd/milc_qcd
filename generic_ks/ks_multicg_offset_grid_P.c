/******* ks_multicg_offset_grid_P.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the multishift GRID inverter */

/* C. DeTar 7/10/17 created */

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

extern GRID_4Dgrid *grid_full;
extern GRID_4DRBgrid *grid_rb;

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
		    GRID_resid_arg_t **grid_resid_arg, int nmass)
{
  /* For now we don't support separate residuals for each mass */
  for(int i=0; i<nmass; ++i){
    qic[i].final_rsq     = grid_resid_arg[i]->final_rsq;
    qic[i].final_relrsq  = 0.;                            /* Not supported at the moment */
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
  if(initialize_grid() != GRID_SUCCESS){
    node0_printf("%s: Error initializing GRID\n",myname);
    terminate(1);
  }

  double tot_cg_time = -dclock();
  double nflop = 1205 + 15*nmass;

  assert(qic[0].parity != EVENANDODD && "EVENANDODD not yet implemented");
    
  if(qic[0].relresid != 0.){
    node0_printf("%s: Grid code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
  }

  if( nmass==0 )return 0;

  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic+0, nmass );

  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic+0, nmass );

  /* Map the masses */
  for(i = 0; i < nmass; i++)
    mass[i] = sqrt(ksp[i].offset/4.0);

  /* Map the input and output fields */
  double t_sp1 = -dclock();
  grid_src = GRID_create_V_from_vec( src, qic[0].parity, grid_full, grid_rb);
  t_sp1 += dclock();
  
  /* Create the solution fields, but leave them zeroed out */
  double t_sp2 = -dclock();
  for(i = 0; i < nmass; i++){
    //    grid_sol[i] =  GRID_create_V( qic[0].parity);
    grid_sol[i] =  GRID_create_V_from_vec(psim[i], qic[0].parity, grid_full, grid_rb);
  }
  t_sp2 += dclock();

  fflush(stdout);

  double t_l = -dclock();
  links = GRID_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full);
  t_l += dclock();
  
  double dtimegridinv = -dclock();
#if 0
  node0_printf("Faking GRID_ks_multicg_offset\n");fflush(stdout);

  for(int i = 0; i < nmass; i++){
    GRID_asqtad_invert( &info, links, &grid_invert_arg, grid_resid_arg[i], 
			mass[i], grid_sol[i], grid_src, grid_full, grid_rb );
    num_iters += grid_resid_arg[i]->final_iter;
  }
#else
  GRID_asqtad_invert_multi( &info, links, &grid_invert_arg, grid_resid_arg,
			    mass, nmass, grid_sol, grid_src, grid_full, grid_rb );
  /* Take the maximum number of iters in order to compare with other multi-mass inverters */
  num_iters = 0;
  for(int i = 0; i < nmass; i++){
    if(grid_resid_arg[i]->final_iter > num_iters)
      num_iters = grid_resid_arg[i]->final_iter;
  }
#endif
  dtimegridinv += dclock();
  double dtimeinv = info.final_sec;
  double t_gr = dtimegridinv - dtimeinv;
  get_grid_resid_arg( qic, grid_resid_arg, nmass);

  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg, nmass);
  
  /* Unpack the solutions */
  double t_sl = -dclock();
  for(i=0; i<nmass; ++i)
    /* Copy results back to su3_vector */
    GRID_extract_V_to_vec( psim[i], grid_sol[i], qic->parity);
  t_sl += dclock();
  
  /* Free GRID fields  */
  
  GRID_destroy_V(grid_src);    
  for(i = 0; i < nmass; i++)
    GRID_destroy_V(grid_sol[i]);     

  GRID_asqtad_destroy_L(links);
  
  tot_cg_time += dclock();

#ifdef CGTIME
  const char *prec_label[2] = {"F", "D"};
  if(this_node==0){
    printf("CONGRAD5: time = %e "
           "(multicg_offset_Grid %s) "
           "masses = %d iters = %d "
           "mflops = %e"
	   "\n",
           tot_cg_time,
           prec_label[qic[0].prec-1],nmass,num_iters,
           (double)(nflop)*volume*num_iters/(1.0e6*tot_cg_time*numnodes())
           );
    fflush(stdout);}
#ifdef REMAP
  if(this_node==0){
    printf("MILC<-->Grid data layout conversion timings\n"
	   "\t src-spinor   = %e\n"
	   "\t dest-spinors = %e\n"
	   "\t soln-spinor  = %e\n"
	   "\t links        = %e\n"
	   "\t Grid remap   = %e\n"
	   "\t ---------------------------\n"
	   "\t total remap  = %e\n"
	   , t_sp1, t_sp2, t_l, t_sl, t_gr
	   , t_sp1 + t_sp2 + t_l + t_sl + t_gr
	   );
    printf("CONGRAD5-Grid: "
  	   "solve-time = %e "
           "(multicg_offset_Grid %s) "
           "masses = %d iters = %d "
           "mflops(ignore data-conv.) = %e\n",
           dtimeinv,
           prec_label[qic[0].prec-1],nmass,num_iters,
           (double)(nflop)*volume*num_iters/(1.0e6*dtimeinv*numnodes())
           );
    fflush(stdout);}
#endif
#endif

  total_iters += num_iters;
  return num_iters;
}
