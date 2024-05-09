/******************* d_congrad5_fn_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* This is the MILC wrapper for the Grid single-mass inverter */

/* NOTE: This code is actually an include file for d_congrad5_fn_grid_F.c
   and d_congrad5_fn_grid_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_CONGRAD_PARITY_GRID
   KS_CONGRAD_MIXED_PARITY_GRID
   KS_CONGRAD_BLOCK_PARITY_GRID
   KS_CONGRAD_MIXED_BLOCK_PARITY_GRID

*/

/* Redefinitions according to requested precision */

#if ( GRID_PrecisionInt == 1 )

#define KS_CONGRAD_PARITY_GRID         ks_congrad_parity_grid_F
#define KS_CONGRAD_BLOCK_PARITY_GRID   ks_congrad_block_parity_grid_F
#define MYREAL GRID_F_Real

#else

#define KS_CONGRAD_PARITY_GRID              ks_congrad_parity_grid_D
#define KS_CONGRAD_BLOCK_PARITY_GRID        ks_congrad_block_parity_grid_D
#define KS_CONGRAD_MIXED_PARITY_GRID        ks_congrad_mixed_parity_grid_D
#define KS_CONGRAD_MIXED_BLOCK_PARITY_GRID  ks_congrad_mixed_block_parity_grid_D
#define MYREAL GRID_D_Real

#endif


#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
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
		     quark_invert_control *qic )
{
  grid_invert_arg->parity       = milc2grid_parity(qic->parity);
  grid_invert_arg->max          = qic->max;
  grid_invert_arg->maxInner     = qic->max_inner;
  grid_invert_arg->nrestart     = qic->nrestart;
}


static GRID_resid_arg_t *
create_grid_resid_arg( quark_invert_control *qic )
{
  GRID_resid_arg_t *res_arg;
  char myname[] = "create_grid_resid_arg";

  /* Pointers for residual errors */
  res_arg = (GRID_resid_arg_t *)malloc(sizeof(GRID_resid_arg_t ));
  if(res_arg == NULL){
    printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
    terminate(1);
  }
  *res_arg = GRID_RESID_ARG_DEFAULT;
  /* For now the residuals are the same for all sources and masses */
  res_arg->resid        = qic->resid;
  res_arg->relresid     = 0.;  /* NOT SUPPORTED */
  res_arg->final_rsq    = 0.;
  res_arg->final_rel    = 0.;

  return res_arg;
}

/* Collect inversion statistics */

static void 
get_grid_resid_arg( quark_invert_control *qic, 
		    GRID_resid_arg_t* grid_resid_arg )
{
  qic->final_rsq     = grid_resid_arg->final_rsq;
  qic->final_relrsq  = 0.;                          /* Not supported at the moment */
  qic->size_r        = grid_resid_arg->size_r;
  qic->size_relr     = grid_resid_arg->size_relr;
  qic->final_iters   = grid_resid_arg->final_iter;
  qic->final_restart = grid_resid_arg->final_restart;
}

static void
destroy_grid_resid_arg(GRID_resid_arg_t *res_arg)
{
  free(res_arg);
}


/*! \brief call to ks_congrad_parity_grid.
 *
 */
int
KS_CONGRAD_PARITY_GRID ( su3_vector *src,
			 su3_vector *sol,
			 quark_invert_control *qic,
			 Real mass,
			 fn_links_t *fn)			 
{
  char myname[] = "ks_congrad_parity_grid";
  double nflop = 1187;
  int iters = 0;
  GRID_info_t info = GRID_INFO_ZERO;
  GRID_ColorVector *grid_sol, *grid_src;
  GRID_resid_arg_t *grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_FermionLinksAsqtad  *links;    
  
  double tot_cg_time = -dclock();
  
  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic );
  
  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic );
  
  /* Data layout conversions */
  
  double t_sp1, t_sp2, t_l;
  t_sp1 = -dclock();
  
  grid_src = GRID_create_V_from_vec( src, qic->parity, grid_full, grid_rb);
  
  t_sp1 += dclock();
  t_sp2 = -dclock();
  
  grid_sol = GRID_create_V_from_vec( sol, qic->parity, grid_full, grid_rb);
  
  t_sp2 += dclock();
  t_l   = -dclock(); 
  
  /* For now we are taking the thin links from the site structure, so the first parameter is NULL */
  links = GRID_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full );

  t_l   += dclock(); 
  
  
  double dtimegridinv = -dclock();
  
  GRID_asqtad_invert( &info, links, &grid_invert_arg, 
		      grid_resid_arg, (MYREAL)mass, grid_sol, grid_src, grid_full, grid_rb );
  iters = grid_resid_arg->final_iter;
  
  dtimegridinv += dclock();
  double dtimeinv = info.final_sec;
  double t_gr = dtimegridinv - dtimeinv - info.misc_sec;
  
  get_grid_resid_arg(qic, grid_resid_arg);
  
  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg);
  
  /* Copy results back to su3_vector */
  double t_sl = -dclock();
  GRID_extract_V_to_vec( sol, grid_sol, qic->parity);
  t_sl += dclock();

#if 0
  /* Fix normalization */
  int i;
  Real rescale = 1./(2.*mass);
  FORSOMEFIELDPARITY_OMP(i, qic->parity, default(shared)){
    scalar_mult_su3_vector( sol+i, rescale, sol+i);
  } END_LOOP_OMP;
#endif
  
  /* Free GRID fields  */
  
  GRID_destroy_V(grid_src);    
  GRID_destroy_V(grid_sol);     
  GRID_asqtad_destroy_L(links);
  
  tot_cg_time +=dclock();

#ifdef CGTIME
  char *prec_label[2] = {"F", "D"};
  if(this_node==0) {
    printf("CONGRAD5: time = %e "
	   "(Grid %s) masses = 1 iters = %d "
	   "mflops = %e "
	   "\n"
	   , tot_cg_time
	   , prec_label[GRID_PrecisionInt-1], iters
	   , (double)(nflop*volume*iters/(1.0e6*tot_cg_time*numnodes()))
	   );
    fflush(stdout);
  }
#ifdef REMAP
  if(this_node==0) {
    printf("MILC<-->Grid data layout conversion timings\n"
	   "\t src-spinor    = %e\n"
	   "\t dest-spinors  = %e\n"
	   "\t soln-spinor   = %e\n"
	   "\t links         = %e\n"
	   "\t Grid overhead = %e\n"
	   "\t ---------------------------\n"
	   "\t total remap   = %e\n"
	   , t_sp1, t_sp2, t_l, t_sl, t_gr
	   , t_sp1 + t_sp2 + t_l + t_sl + t_gr
	   );
    printf("CONGRAD5-Grid: "
	   "solve-time = %e "
	   "(Grid %s) masses = 1 iters = %d "
	   "mflops(ignore data-conv.) = %e "
	   "\n"
	   , dtimeinv
	   , prec_label[GRID_PrecisionInt-1], iters
	   , (double)(nflop*volume*iters/(1.0e6*dtimeinv*numnodes()))
	   );
    fflush(stdout);
  }
#endif
#endif

  return iters;
}


/*! \brief call to ks_congrad_mixed_parity_grid.
 *
 */
int
KS_CONGRAD_MIXED_PARITY_GRID ( su3_vector *src,
			       su3_vector *sol,
			       quark_invert_control *qic,
			       Real mass,
			       fn_links_t *fn)			 
{
  char myname[] = "ks_congrad_mixed_parity_grid";
  double nflop = 1187;
  int iters = 0;
  GRID_info_t info = GRID_INFO_ZERO;
  GRID_D3_ColorVector *grid_sol, *grid_src;
  GRID_resid_arg_t *grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_D3_FermionLinksAsqtad  *links_d;    
  GRID_F3_FermionLinksAsqtad  *links_f;    
  
  double tot_cg_time = -dclock();
  
  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic );
  
  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic );
  
  /* Data layout conversions */
  
  double t_sp1, t_sp2, t_l;
  t_sp1 = -dclock();
  
  grid_src = GRID_D3_create_V_from_vec( src, qic->parity, grid_full, grid_rb);
  
  t_sp1 += dclock();
  t_sp2 = -dclock();
  
  grid_sol = GRID_D3_create_V_from_vec( sol, qic->parity, grid_full, grid_rb);
  
  t_sp2 += dclock();
  t_l   = -dclock(); 
  /* For now we are taking the thin links from the site structure, so the first parameter is NULL */
  links_d = GRID_D3_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full );
  links_f = GRID_F3_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full );
  t_l   += dclock(); 
  
  double dtimegridinv = -dclock();
  GRID_D3_asqtad_invert_mixed( &info, links_f, links_d, &grid_invert_arg, 
			       grid_resid_arg, (double)mass, grid_sol, grid_src, grid_full, grid_rb );
  iters = grid_resid_arg->final_iter;
  
  dtimegridinv += dclock();
  double dtimeinv = info.final_sec;
  double t_gr = dtimegridinv - dtimeinv - info.misc_sec;
  
  get_grid_resid_arg(qic, grid_resid_arg);
  
  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg);
  
  /* Copy results back to su3_vector */
  double t_sl = -dclock();
  GRID_D3_extract_V_to_vec( sol, grid_sol, qic->parity);
  t_sl += dclock();

#if 0
  /* Fix normalization */
  int i;
  Real rescale = 1./(2.*mass);
  FORSOMEFIELDPARITY_OMP(i, qic->parity, default(shared)){
    scalar_mult_su3_vector( sol+i, rescale, sol+i);
  } END_LOOP_OMP;
#endif
  
  /* Free GRID fields  */
  
  GRID_D3_destroy_V(grid_src);    
  GRID_D3_destroy_V(grid_sol);     
  GRID_D3_asqtad_destroy_L(links_d);
  GRID_F3_asqtad_destroy_L(links_f);
  
  tot_cg_time +=dclock();

#ifdef CGTIME
  char *prec_label[2] = {"F", "D"};
  if(this_node==0) {
    printf("CONGRAD5: time = %e "
	   "(Grid %s) masses = 1 iters = %d "
	   "mflops = %e "
	   "\n"
	   , tot_cg_time
	   , prec_label[GRID_PrecisionInt-1], iters
	   , (double)(nflop*volume*iters/(1.0e6*tot_cg_time*numnodes()))
	   );
    fflush(stdout);
  }
#ifdef REMAP
  if(this_node==0) {
    printf("MILC<-->Grid data layout conversion timings\n"
	   "\t src-spinor    = %e\n"
	   "\t dest-spinors  = %e\n"
	   "\t soln-spinor   = %e\n"
	   "\t links         = %e\n"
	   "\t Grid overhead = %e\n"
	   "\t ---------------------------\n"
	   "\t total remap   = %e\n"
	   , t_sp1, t_sp2, t_l, t_sl, t_gr
	   , t_sp1 + t_sp2 + t_l + t_sl + t_gr
	   );
    printf("CONGRAD5-Grid: "
	   "solve-time = %e "
	   "(Grid D mixed) masses = 1 iters = %d "
	   "mflops(ignore data-conv.) = %e "
	   "\n"
	   , dtimeinv
	   , iters
	   , (double)(nflop*volume*iters/(1.0e6*dtimeinv*numnodes()))
	   );
    fflush(stdout);
  }
#endif
#endif

  return iters;
}


/*! \brief call to ks_congrad_block_parity_grid.
 *
 */
int
KS_CONGRAD_BLOCK_PARITY_GRID ( int nrhs,
			       su3_vector *src[],
			       su3_vector *sol[],
			       quark_invert_control *qic,
			       Real mass,
			       fn_links_t *fn)			 
{
  char myname[] = "ks_congrad_block_parity_grid";
  double nflop = 1187;
  int iters = 0;
  GRID_info_t info = GRID_INFO_ZERO;
  GRID_ColorVectorBlock *grid_sol, *grid_src;
  GRID_resid_arg_t *grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_FermionLinksAsqtad  *links;    
  
  double tot_cg_time = -dclock();
  
  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic );
  
  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic );
  
  /* Data layout conversions */

  GRID_5Dgrid *grid_5D = GRID_create_5Dgrid(nrhs, grid_full);
  GRID_5DRBgrid *grid_5Drb = GRID_create_5DRBgrid(nrhs, grid_full);
  
  double t_sp1, t_sp2, t_l;
  t_sp1 = -dclock();
  
  grid_src = GRID_create_nV_from_vecs( src, nrhs, qic->parity, 
				       grid_5D, grid_5Drb, grid_full, grid_rb);
  
  t_sp1 += dclock();
  t_sp2 = -dclock();
  
  grid_sol = GRID_create_nV_from_vecs( sol, nrhs, qic->parity,
				       grid_5D, grid_5Drb, grid_full, grid_rb);
  
  t_sp2 += dclock();
  t_l   = -dclock(); 
  
  /* For now we are taking the thin links from the site structure, so the first parameter is NULL */
  links = GRID_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full);

  t_l   += dclock(); 
  
  
  double dtimegridinv = -dclock();
  
  GRID_asqtad_invert_block( &info, links, &grid_invert_arg, 
			    grid_resid_arg, (MYREAL)mass, nrhs, grid_sol, grid_src,
			    grid_5D, grid_5Drb, grid_full, grid_rb);

  iters = grid_resid_arg->final_iter;
  
  dtimegridinv += dclock();
  double dtimeinv = info.final_sec;
  double t_gr = dtimegridinv - dtimeinv - info.misc_sec;
  
  get_grid_resid_arg(qic, grid_resid_arg);
  
  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg);
  
  /* Copy results back to su3_vector */
  double t_sl = -dclock();
  GRID_extract_nV_to_vecs( sol, nrhs, grid_sol, qic->parity);
  t_sl += dclock();

#if 0
  /* Fix normalization */
  int i;
  Real rescale = 1./(2.*mass);
  FORSOMEFIELDPARITY_OMP(i, qic->parity, default(shared)){
    for(int j = 0; j < nrhs; j++)
      scalar_mult_su3_vector( sol[j]+i, rescale, sol[j]+i);
  } END_LOOP_OMP;
#endif
  
  /* Free GRID fields  */
  
  GRID_destroy_nV(grid_src);    
  GRID_destroy_nV(grid_sol);     
  GRID_asqtad_destroy_L(links);
  
  /* Free 5D lattice */
  GRID_destroy_5Dgrid(grid_5D);
  GRID_destroy_5DRBgrid(grid_5Drb);

  tot_cg_time +=dclock();

#ifdef CGTIME
  char *prec_label[2] = {"F", "D"};
  if(this_node==0) {
    printf("CONGRAD5: time = %e "
	   "(Grid-block %s) masses = 1 iters = %d rhs = %d "
	   "mflops = %e "
	   "\n"
	   , tot_cg_time
	   , prec_label[GRID_PrecisionInt-1], iters, nrhs
	   , (double)(nflop*nrhs*volume*iters/(1.0e6*tot_cg_time*numnodes()))
	   );
    fflush(stdout);
  }
#ifdef REMAP
  if(this_node==0) {
    printf("MILC<-->Grid data layout conversion timings\n"
	   "\t src-spinor    = %e\n"
	   "\t dest-spinors  = %e\n"
	   "\t soln-spinor   = %e\n"
	   "\t links         = %e\n"
	   "\t Grid overhead = %e\n"
	   "\t ---------------------------\n"
	   "\t total remap   = %e\n"
	   , t_sp1, t_sp2, t_l, t_sl, t_gr
	   , t_sp1 + t_sp2 + t_l + t_sl + t_gr
	   );
    printf("CONGRAD5-Grid: "
	   "solve-time = %e "
	   "(Grid %s) masses = 1 iters = %d rhs = %d "
	   "mflops(ignore data-conv.) = %e "
	   "\n"
	   , dtimeinv
	   , prec_label[GRID_PrecisionInt-1], iters, nrhs
	   , (double)(nflop*nrhs*volume*iters/(1.0e6*dtimeinv*numnodes()))
	   );
    fflush(stdout);
  }
#endif
#endif

  return iters;
}

/*! \brief call to ks_congrad_mixed_block_parity_grid.
 *
 */
int
KS_CONGRAD_MIXED_BLOCK_PARITY_GRID ( int nrhs,
				     su3_vector *src[],
				     su3_vector *sol[],
				     quark_invert_control *qic,
				     Real mass,
				     fn_links_t *fn)			 
{
  char myname[] = "ks_congrad_mixed_block_parity_grid";
  double nflop = 1187;
  int iters = 0;
  GRID_info_t info = GRID_INFO_ZERO;
  GRID_D3_ColorVectorBlock *grid_sol, *grid_src;
  GRID_resid_arg_t *grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_D3_FermionLinksAsqtad  *links_d;    
  GRID_F3_FermionLinksAsqtad  *links_f;    
  
  double tot_cg_time = -dclock();
  
  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic );
  
  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic );
  
  /* Data layout conversions */

  GRID_5Dgrid *grid_5D = GRID_create_5Dgrid(nrhs, grid_full);
  GRID_5DRBgrid *grid_5Drb = GRID_create_5DRBgrid(nrhs, grid_full);

  double t_sp1, t_sp2, t_l;
  t_sp1 = -dclock();
  
  grid_src = GRID_D3_create_nV_from_vecs( src, nrhs, qic->parity, 
					  grid_5D, grid_5Drb, grid_full, grid_rb);
  
  t_sp1 += dclock();
  t_sp2 = -dclock();
  
  grid_sol = GRID_D3_create_nV_from_vecs( sol, nrhs, qic->parity,
					  grid_5D, grid_5Drb, grid_full, grid_rb);
  
  t_sp2 += dclock();
  t_l   = -dclock(); 

  /* For now we are taking the thin links from the site structure, so the first parameter is NULL */
  links_d = GRID_D3_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full);
  links_f = GRID_F3_asqtad_create_L_from_MILC( NULL, get_fatlinks(fn), get_lnglinks(fn), grid_full);

  t_l   += dclock(); 
  
  
  double dtimegridinv = -dclock();
  
  GRID_D3_asqtad_invert_mixed_block( &info, links_f, links_d, &grid_invert_arg, 
				     grid_resid_arg, (double)mass, nrhs,
				     grid_sol, grid_src,
				     grid_5D, grid_5Drb, grid_full, grid_rb);
  
  iters = grid_resid_arg->final_iter;
  
  dtimegridinv += dclock();
  double dtimeinv = info.final_sec;
  double t_gr = dtimegridinv - dtimeinv - info.misc_sec;
  
  get_grid_resid_arg(qic, grid_resid_arg);
  
  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg);
  
  /* Copy results back to su3_vector */
  double t_sl = -dclock();
  GRID_D3_extract_nV_to_vecs( sol, nrhs, grid_sol, qic->parity);
  t_sl += dclock();

#if 0
  /* Fix normalization */
  int i;
  Real rescale = 1./(2.*mass);
  FORSOMEFIELDPARITY_OMP(i, qic->parity, default(shared)){
    for(int j = 0; j < nrhs; j++)
      scalar_mult_su3_vector( sol[j]+i, rescale, sol[j]+i);
  } END_LOOP_OMP;
#endif
  
  /* Free GRID fields  */
  
  GRID_D3_destroy_nV(grid_src);    
  GRID_D3_destroy_nV(grid_sol);     
  GRID_D3_asqtad_destroy_L(links_d);
  GRID_F3_asqtad_destroy_L(links_f);
  
  /* Free 5D lattice */
  GRID_destroy_5Dgrid(grid_5D);
  GRID_destroy_5DRBgrid(grid_5Drb);

  tot_cg_time +=dclock();

#ifdef CGTIME
  char *prec_label[2] = {"F", "D"};
  if(this_node==0) {
    printf("CONGRAD5: time = %e "
	   "(Grid-block %s) masses = 1 iters = %d rhs = %d "
	   "mflops = %e "
	   "\n"
	   , tot_cg_time
	   , prec_label[GRID_PrecisionInt-1], iters, nrhs
	   , (double)(nflop*nrhs*volume*iters/(1.0e6*tot_cg_time*numnodes()))
	   );
    fflush(stdout);
  }
#ifdef REMAP
  if(this_node==0) {
    printf("MILC<-->Grid data layout conversion timings\n"
	   "\t src-spinor    = %e\n"
	   "\t dest-spinors  = %e\n"
	   "\t soln-spinor   = %e\n"
	   "\t links         = %e\n"
	   "\t Grid overhead = %e\n"
	   "\t ---------------------------\n"
	   "\t total remap   = %e\n"
	   , t_sp1, t_sp2, t_l, t_sl, t_gr
	   , t_sp1 + t_sp2 + t_l + t_sl + t_gr
	   );
    printf("CONGRAD5-Grid: "
	   "solve-time = %e "
	   "(Grid %s) masses = 1 iters = %d rhs = %d "
	   "mflops(ignore data-conv.) = %e "
	   "\n"
	   , dtimeinv
	   , prec_label[GRID_PrecisionInt-1], iters, nrhs
	   , (double)(nflop*nrhs*volume*iters/(1.0e6*dtimeinv*numnodes()))
	   );
    fflush(stdout);
  }
#endif
#endif

  return iters;
}


