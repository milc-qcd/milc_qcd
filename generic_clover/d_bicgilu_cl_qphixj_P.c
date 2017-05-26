/******* d_bicgilu_cl_qphixj_P.c - BiCG ILU ***********************/
/* MIMD version 7 */
/* 4/27/17 C. DeTar */
/* Modified to map straight from MILC to QPhiX without going through raw fields */

/* This is the MILC wrapper for the Balint Joo's QPhiX clover inverter */
/* MILC fields are in the MILC generic precision.  They are converted to
   the desired QPHIX precision for the inner solve */

/* NOTE: This code is actually an include file for d_bicgilu_cl_qphixj_F.c
   and d_bicgilu_cl_qphixj_D.c, so any edits should be consistent with this
   purpose. 

*/    

#include "generic_clover_includes.h"
#include "../include/generic_qphixj.h"
#include "../include/openmp_defs.h"
#include "../include/qphixj/qphixj.h"

/* Redefinitions according to requested precision */

#if ( QPHIXJ_PrecisionInt == 1 )
#define BICGILU_CL_QPHIXJ_INNER bicgilu_cl_qphixj_inner_F
#else
#define BICGILU_CL_QPHIXJ_INNER bicgilu_cl_qphixj_inner_D
#endif

/* adjoint of an SU3 matrix with no precision conversion */
/* result in b */
static void 
su3_adjoint_p( su3_matrix *a, su3_matrix *b )
{
  register int i,j;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
  CONJG( a->e[j][i], b->e[i][j] );
  }
}

/*!
 * Copy backlinks with adjoint
 */
static su3_matrix *
create_backlinks_with_adjoint_from_site(void)
{
  su3_matrix *t_bl = NULL;
  register int i;
  register site *s;
  int dir;
  msg_tag *tag[4];
  char myname[] = "create_backlinks_with_adjoint_from_site";
  
  /* Allocate space for t_lbl */
  t_bl = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_bl==NULL){
    printf("%s(%d): no room for t_lbl\n",myname,this_node);
    terminate(1);
  }
  
  /* gather backwards gauge links */
  for( dir=XUP; dir<=TUP; dir ++){
    tag[dir] = start_gather_site( F_OFFSET(link[dir])
				  , sizeof(su3_matrix)
				  , OPP_DIR(dir)
				  , EVENANDODD
				  , gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLFIELDSITES_OMP(i,) {
      su3_adjoint_p((su3_matrix *)gen_pt[dir][i], t_bl + dir + 4*i);
    } END_LOOP_OMP;
    cleanup_gather( tag[dir] );
  }
  
  return t_bl;
}

static void
destroy_backlinks(su3_matrix *t_bl)
{
  free(t_bl);
}

void BICGILU_CL_QPHIXJ_INNER(clover *milc_clov, Real kappa, wilson_vector r[], 
			     wilson_vector dest[], quark_invert_control *qic)
{

  /* Map clover term to raw with possible precision conversion */
  char myname[] = "bicgilu_cl_qphixj_inner";
  double dtime = -dclock();

  initialize_qphixj();

  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time to initialize qphixj %.2e\n", dtime);
#endif
  
  dtime = -dclock();
  /* Backward gauge links with possible precision conversion */
  su3_matrix* fieldback = create_backlinks_with_adjoint_from_site();

  /* Create QPHIX objects */
  QPHIXJ_FermionLinksWilson  *wilson = 
    QPHIXJ_wilson_create_L_from_MILC( fieldback, milc_clov, kappa, QPHIXJ_EVEN );

  destroy_backlinks(fieldback);

  QPHIXJ_DiracFermion *in = QPHIXJ_create_D_from_wvec( r, QPHIXJ_EVEN );
  QPHIXJ_DiracFermion *out = QPHIXJ_create_D_from_wvec( dest, QPHIXJ_EVEN );
  
  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time to create back links and import fields %.2e\n", dtime);
#endif
  QPHIXJ_invert_arg_t inv_arg;
  inv_arg.max  = qic->max*qic->nrestart;
  inv_arg.restart = qic->max;
  inv_arg.nrestart = qic->nrestart;
  
  QPHIXJ_info_t info = QPHIXJ_INFO_ZERO;
  
  QPHIXJ_resid_arg_t res_arg = QPHIXJ_RESID_ARG_DEFAULT;
  res_arg.resid = qic->resid*qic->resid;
  res_arg.relresid = qic->relresid*qic->relresid;

  /* Do the inversion on the even sites. Note: kappa is ignored here.  See create_L above */
  QPHIXJ_wilson_invert( &info, wilson, &inv_arg, &res_arg,
			kappa, out, in);

  dtime = -dclock();
  /* Map QPhiXJ out to MILC dest with possible precision conversion */
  QPHIXJ_extract_D_to_wvec( dest, out, EVEN );
  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time to export fields %.2e\n", dtime);
#endif

  /* Convergence data */
  qic->final_rsq     = res_arg.final_rsq;
  qic->final_relrsq  = res_arg.final_rel;
  qic->final_iters   = res_arg.final_iter;
  qic->final_restart = res_arg.final_restart;
  qic->converged = (res_arg.final_rsq <= qic->resid) ? 1 : 0;
  
  /* Relative residual. Not supported here */
  qic->size_r = 0.0;
  qic->size_relr = 0.0;
  
  /* Clean up */
  
  QPHIXJ_wilson_destroy_L(wilson);
  QPHIXJ_destroy_D( out );
  QPHIXJ_destroy_D( in );

}
