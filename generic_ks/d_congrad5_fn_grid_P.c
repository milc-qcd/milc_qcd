/******************* d_congrad5_fn_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* 11/28/15 Created by Dhiraj Khalamkar */
/* 12/23/15 Modified by C. DeTar */

/* This is the MILC wrapper for the Grid single-mass inverter */

/* NOTE: This code is actually an include file for d_congrad5_fn_grid_F.c
   and d_congrad5_fn_grid_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_CONGRAD_PARITY_GRID

*/

/* Redefinitions according to requested precision */

#if ( GRID_PrecisionInt == 1 )

#define KS_CONGRAD_PARITY_GRID   ks_congrad_parity_grid_F
#define MYREAL GRID_F_Real
#define MYSU3_MATRIX fsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_F_G

#else

#define KS_CONGRAD_PARITY_GRID   ks_congrad_parity_grid_D
#define MYREAL GRID_D_Real
#define MYSU3_MATRIX dsu3_matrix
#define COPY_MILC_TO_G copy_milc_to_D_G

#endif

#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <assert.h>

/* Conversions from prevailing MILC formats to specified formats */

#if (PRECISION==1)

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

#else

/* Convert (or copy) su3_matrix from prevailing to single precision */
static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

#endif

#define copy_milc_to_F_G(d,s) p2f_mat(d,s);
#define copy_milc_to_D_G(d,s) p2d_mat(d,s);

/* Load inversion args for Level 3 inverter */

static void 
set_grid_invert_arg( GRID_invert_arg_t* grid_invert_arg, 
		      quark_invert_control *qic )
{
  grid_invert_arg->parity       = milc2grid_parity(qic->parity);
  grid_invert_arg->max          = qic->max;
  grid_invert_arg->nrestart     = qic->nrestart;
}


static GRID_resid_arg_t *
create_grid_resid_arg( quark_invert_control *qic )
{
  GRID_resid_arg_t *res_arg;
  char myname[] = "create_grid_resid_arg";
  int isrc,imass;

  /* Pointers for residual errors */
  res_arg = (GRID_resid_arg_t *)malloc(sizeof(GRID_resid_arg_t ));
  if(res_arg == NULL){
    printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
    terminate(1);
  }
  *res_arg = GRID_RESID_ARG_DEFAULT;
  /* For now the residuals are the same for all sources and masses */
  res_arg->resid = qic->resid * qic->resid;
  res_arg->relresid     = 0.;  /* NOT SUPPORTED */
  res_arg->final_rsq    = 0.;
  res_arg->final_rel    = 0.;

  return res_arg;
}

/* Collect inversion statistics */

static void 
get_grid_resid_arg( quark_invert_control *qic, 
		     GRID_resid_arg_t* grid_resid_arg, int iters )
{
  qic->final_rsq     = grid_resid_arg->final_rsq;
  qic->final_relrsq  = 0.;                          /* Not supported at the moment */
  qic->size_r        = grid_resid_arg->size_r;
  qic->size_relr     = grid_resid_arg->size_relr;
  qic->final_iters   = iters;
  qic->final_restart = grid_resid_arg->final_restart;
}

static void
destroy_grid_resid_arg(GRID_resid_arg_t *res_arg)
{
  free(res_arg);
}


/*!
 * Copy fatbacklinks without the adjoint.
 * Grid does the adjoint in the generated code for the dslash kernels for the back
 * links.
 */
static MYSU3_MATRIX *
create_fatbacklinks_without_adjoint(su3_matrix *t)
{
  MYSU3_MATRIX *t_bl = NULL;
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "create_fatbacklinks_without_adjoint";
  
  /* Allocate space for t_lbl */
  t_bl = (MYSU3_MATRIX *)malloc(sites_on_node*4*sizeof(MYSU3_MATRIX));
  if(t_bl==NULL){
    printf("%s(%d): no room for t_lbl\n",myname,this_node);
    terminate(1);
  }
  
  tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
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
      COPY_MILC_TO_G(t_bl + dir + 4*i, (su3_matrix *)gen_pt[dir][i]);
    } END_LOOP_OMP
    cleanup_gather( tag[dir] );
  }
  
  free(tempmat1); 
  tempmat1 = NULL;
  return t_bl;
}

/*!
 * Copy lngbacklinks without the adjoint.
 * Grid does the adjoint in the generated code for the dslash kernels for the back
 * links.
 */
static MYSU3_MATRIX *
create_lngbacklinks_without_adjoint(su3_matrix *t)
{
  MYSU3_MATRIX *t_bl = NULL;
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_backlinks_without_adjoint";
  
  /* Allocate space for t_lbl if NULL */
  t_bl = (MYSU3_MATRIX *)malloc(sites_on_node*4*sizeof(MYSU3_MATRIX));
  if(t_bl==NULL){
    printf("%s(%d): no room for t_lbl\n",myname,this_node);
    terminate(1);
  }
  
  tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }
  
  /* gather backwards longlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLFIELDSITES_OMP(i,){
      tempmat1[i] = t[dir+4*i];
    } END_LOOP_OMP
    tag[dir] = start_gather_field( tempmat1
				   , sizeof(su3_matrix)
				   , OPP_3_DIR(DIR3(dir))
				   , EVENANDODD
				   , gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLFIELDSITES_OMP(i,) {
      MYSU3_MATRIX * temp_ = (t_bl + dir + 4*i);
      COPY_MILC_TO_G(t_bl + dir + 4*i, (su3_matrix *)gen_pt[dir][i]);
    } END_LOOP_OMP
    cleanup_gather( tag[dir] );
  }
  
  free(tempmat1); 
  tempmat1 = NULL;
  
  return t_bl;
}

static void
destroy_backlinks(MYSU3_MATRIX *t_bl){
  free(t_bl);
}

/*! 
 * Create the Grid fermion link structure from forward FN links
 */
GRID_FermionLinksAsqtad *
create_grid_L_from_fn_links(fn_links_t *fn, int parity){
  MYSU3_MATRIX *raw_fat_links, *raw_lng_links, *raw_fatback_links, *raw_lngback_links;
  static GRID_FermionLinksAsqtad *links;
  
  raw_fat_links  = create_grid_raw4_G_from_field(get_fatlinks(fn), parity);
  if(raw_fat_links == NULL)terminate(1);

  raw_lng_links = create_grid_raw4_G_from_field(get_lnglinks(fn), parity);
  if(raw_lng_links == NULL)terminate(1);
  
  //  if(get_fatbacklinks(fn) == NULL)
  //  raw_fatback_links = create_fatbacklinks_without_adjoint(get_fatlinks(fn));
  // else
  // WE MUST UNDO THE ADJOINT HERE WHEN WE COPY FROM fn
  // raw_fatback_links = create_grid_raw4_G_from_field(get_fatbacklinks(fn), parity);
  //}
  
  //if(get_lngbacklinks(fn) == NULL)
  // raw_lngback_links = create_lngbacklinks_without_adjoint(get_lnglinks(fn));
  //else
  // WE MUST UNDO THE ADJOINT HERE WHEN WE COPY FROM fn
  //  raw_lngback_links = create_grid_raw4_G_from_field(get_lngbacklinks(fn), parity);
  
  links = GRID_asqtad_create_L_from_raw( (MYREAL *)raw_fat_links, (MYREAL *)raw_lng_links, parity);
  
  destroy_grid_raw4_G(raw_fat_links);
  destroy_grid_raw4_G(raw_lng_links);
  //  destroy_backlinks(raw_fatback_links);
  //  destroy_backlinks(raw_lngback_links);
  
  return links;
}


/*! \brief call to the grid_ks_congrad_parity.
 *
 * Expects that setup_mbench has been called already, so that we do not have to 
 * malloc things agian.
 */
int
KS_CONGRAD_PARITY_GRID ( su3_vector *src
			  , su3_vector *sol
			  , quark_invert_control *qic
			  , Real mass
			  , fn_links_t *fn)			 
{
  char myname[] = "ks_congrad_parity_grid";
  double ttime, dctime, tot_cg_time;
  double dtime;
  double nflop = 1187;
  int iters = 0;
  int otherparity;
  GRID_info_t info = GRID_INFO_ZERO;
  GRID_ColorVector *grid_sol, *grid_src;
  GRID_resid_arg_t  *grid_resid_arg;
  GRID_invert_arg_t grid_invert_arg = GRID_INVERT_ARG_DEFAULT;
  GRID_FermionLinksAsqtad  *links;    
  
  assert(qic->parity != EVENANDODD && "EVENANDODD not yet implemented");
  
#ifdef CGTIME
  tot_cg_time = -dclock();
#endif   
  
  /* Initialize GRID if not already done */
  if(initialize_grid(GRID_PrecisionInt) != GRID_SUCCESS){
    node0_printf("%s: Error initializing GRID\n",myname);
    terminate(1);
  }
#ifdef CG_DEBUG
  dctime = -dclock();
#endif
  
  /* Set grid_invert_arg */
  set_grid_invert_arg( & grid_invert_arg, qic );
  
  /* Pointers for residual errors */
  grid_resid_arg = create_grid_resid_arg( qic );
  
  /* Data layout conversions */
  
#ifdef REMAP
  double t_sp1, t_sp2, t_l;
  t_sp1 = -dclock();
#endif    
  
  grid_src = create_grid_V_from_field( src, qic->parity);
  
#ifdef REMAP
  t_sp1 += dclock();
  t_sp2 = -dclock();
#endif    
  
  grid_sol = create_grid_V_from_field( sol, qic->parity);
  
#ifdef REMAP
  t_sp2 += dclock();
  t_l   = -dclock(); 
#endif     
  
  links = create_grid_L_from_fn_links( fn, EVENANDODD );
  
#if CG_DEBUG
  t_l   += dclock(); 
#endif       
  
#ifdef REMAP
  node0_printf("MILC-->Grid data layout conversion timings"
	       " (Unoptimized Gathers).\n"
	       "\t src-spinor  = %e\n"
	       "\t dest-spinor = %e\n"
	       "\t links       = %e\n"
	       "\t total       = %e\n"
	       , t_sp1, t_sp2, t_l
	       , t_sp1 + t_sp2 + t_l
	       );
  fflush(stdout);
#endif
  
#ifdef CG_DEBUG
  dctime +=dclock();
  dtime = -dclock();
#endif    
  
#ifdef CG_DEBUG
  node0_printf("Calling GRID_asqtad_invert\n");fflush(stdout);
#endif

  iters = GRID_asqtad_invert( &info, links, &grid_invert_arg, 
			      grid_resid_arg, (MYREAL)mass, grid_sol, grid_src );
  
#ifdef CG_DEBUG    
  dtime += dclock();
#endif
  
  get_grid_resid_arg(qic, grid_resid_arg, iters);
  
  /* Free the structure */
  destroy_grid_resid_arg(grid_resid_arg);
  
#ifdef CG_DEBUG
  ttime = -dclock();
#endif
  
  /* Copy results back to su3_vector */
  unload_grid_V_to_field( sol, grid_sol, qic->parity);
  
  /* Free GRID fields  */
  
  GRID_destroy_V(grid_src);    
  GRID_destroy_V(grid_sol);     
  GRID_asqtad_destroy_L(links);
  
#ifdef CG_DEBUG
  ttime +=dclock();
  dctime +=ttime;
#endif
#ifdef CGTIME
  tot_cg_time +=dclock();
  if(this_node==0) {
    char *prec_label[2] = {"F", "D"};
    node0_printf("Grid-CONGRAD5: total cg-time = %e "
		 
#ifdef CG_DEBUG
		 "solve-time = %e "
		 "layout-conversion-time = %e "
#endif            
		 "(fn %s) masses = 1 iters = "
		 "%d mflops = %e "
#ifdef CG_DEBUG
		 "mflops(ignore data-conv.) = %e"
#endif
		 "\n"
		 , tot_cg_time
#ifdef CG_DEBUG
		 , dtime, dctime
#endif
		 , prec_label[PRECISION-1], iters
		 , (double)(nflop*volume*iters/(1.0e6*tot_cg_time*numnodes()))
#ifdef CG_DEBUG
		 , (double)(nflop*volume*iters/(1.0e6*dtime*numnodes()))
#endif
		 );
    fflush(stdout);
  }
#endif

  return iters;
}


