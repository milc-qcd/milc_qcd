/******* d_bicgilu_cl_qphixj_P.c - BiCG ILU ***********************/
/* MIMD version 7 */
/* 4/27/17 C. DeTar */

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

/* Redefinitions according to requested precision */

#if ( QPHIXJ_PrecisionInt == 1 )

#define BICGILU_CL_QPHIXJ_INNER bicgilu_cl_qphixj_inner_F
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix
#define MYWILSON_VECTOR fwilson_vector

#else

#define BICGILU_CL_QPHIXJ_INNER bicgilu_cl_qphixj_inner_D
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix
#define MYWILSON_VECTOR dwilson_vector

#endif

static void  
map_milc_clov_to_qphixj_raw(MYREAL *raw_clov, clover *milc_clov, double scale){ 
  int i,j;								
  MYREAL *r;
  
  FOREVENFIELDSITES_OMP(i, private(r,j)){	
    
    r = raw_clov + 72*i;
    for(j=0; j<6; j++){
      r[j]    = milc_clov->clov_diag[i].di[0][j]; 
      r[j+36] = milc_clov->clov_diag[i].di[1][j]; 
    } 
    for(j=0; j<15; j++){
      r[2*j+6]  = milc_clov->clov[i].tr[0][j].real; 
      r[2*j+7]  = milc_clov->clov[i].tr[0][j].imag; 
      r[2*j+42] = milc_clov->clov[i].tr[1][j].real; 
      r[2*j+43] = milc_clov->clov[i].tr[1][j].imag; 
    } 
  } END_LOOP_OMP; 

  FORODDFIELDSITES_OMP(i, private(r,j)){	

    r = raw_clov + 72*i;
    for(j=0; j<6; j++){
      r[j]    = milc_clov->clov_diag[i].di[0][j]*scale; 
      r[j+36] = milc_clov->clov_diag[i].di[1][j]*scale; 
    } 
    for(j=0; j<15; j++){
      r[2*j+6]  = milc_clov->clov[i].tr[0][j].real*scale; 
      r[2*j+7]  = milc_clov->clov[i].tr[0][j].imag*scale; 
      r[2*j+42] = milc_clov->clov[i].tr[1][j].real*scale; 
      r[2*j+43] = milc_clov->clov[i].tr[1][j].imag*scale; 
    } 
  } END_LOOP_OMP;
}

/* adjoint of an SU3 matrix with no precision conversion */
/* result in b */
static void 
su3_adjoint_p( su3_matrix *a, MYSU3_MATRIX *b )
{
  register int i,j;
  for(i=0;i<3;i++)for(j=0;j<3;j++){
  CONJG( a->e[j][i], b->e[i][j] );
  }
}


/*!
 * Copy backlinks with adjoint
 */
static MYSU3_MATRIX *
create_backlinks_with_adjoint_from_site(void)
{
  MYSU3_MATRIX *t_bl = NULL;
  register int i;
  register site *s;
  int dir;
  msg_tag *tag[4];
  char myname[] = "create_backlinks_with_adjoint_from_site";
  
  /* Allocate space for t_lbl */
  t_bl = (MYSU3_MATRIX *)malloc(sites_on_node*4*sizeof(MYSU3_MATRIX));
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
destroy_backlinks(MYSU3_MATRIX *t_bl)
{
  free(t_bl);
}

void BICGILU_CL_QPHIXJ_INNER(clover *milc_clov, Real kappa, wilson_vector r[], 
			     wilson_vector dest[], quark_invert_control *qic)
{

  /* Map clover term to raw with possible precision conversion */
  char myname[] = "bicgilu_cl_qphixj_inner";
  double dtime = -dclock();

  MYREAL *raw_clov     = (MYREAL*)malloc(72*sites_on_node*sizeof(MYREAL));
  if(raw_clov == NULL){
    printf("%s(%d): no room for raw_clov\n", myname, this_node);
    terminate(1);
  }
  double MKsq = -kappa*kappa;
  map_milc_clov_to_qphixj_raw(raw_clov, milc_clov, -4.*MKsq);
  
  initialize_qphixj();

  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time to initialize qphixj %.2e\n", dtime);
#endif
  
  dtime = -dclock();
  /* Backward gauge links with possible precision conversion */
  MYSU3_MATRIX* fieldback = create_backlinks_with_adjoint_from_site();

  /* Create QPHIX objects */
  QPHIXJ_FermionLinksWilson  *wilson = 
    create_qphixj_L_from_fieldback_and_sites( raw_clov, fieldback, QPHIXJ_EVEN );

  destroy_backlinks(fieldback);
  free(raw_clov);
  QPHIXJ_DiracFermion *in = create_qphixj_D_from_field( r, QPHIXJ_EVEN );
  QPHIXJ_DiracFermion *out = create_qphixj_D_from_field( dest, QPHIXJ_EVEN );
  
  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time to import fields %.2e\n", dtime);
#endif
  QPHIXJ_invert_arg_t inv_arg;
  inv_arg.max  = qic->max*qic->nrestart;
  inv_arg.restart = qic->max;
  inv_arg.nrestart = qic->nrestart;
  
  QPHIXJ_info_t info = QPHIXJ_INFO_ZERO;
  
  QPHIXJ_resid_arg_t res_arg = QPHIXJ_RESID_ARG_DEFAULT;
  res_arg.resid = qic->resid*qic->resid;
  res_arg.relresid = qic->relresid*qic->relresid;

  /* Do the inversion on the even sites */
  QPHIXJ_wilson_invert( &info, wilson, &inv_arg, &res_arg,
			kappa, out, in);

  dtime = -dclock();
  /* Map QPhiXJ out to MILC dest with possible precision conversion */
  unload_qphixj_D_to_field( dest, out, EVEN );
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
