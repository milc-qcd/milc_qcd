/******* ks_eigcg_fn_gpu.c - eigCG solver for SU3/fermions ****/
/* MIMD version 7 */

/* GPU version of eigcg.  Can be compiled together with it. */
// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code

#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"
#include <string.h>
#include "../include/generic.h"

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code
int ks_eigcg_parity_gpu(int curr_rhs_idx, int last_rhs_flag, su3_vector *t_src, su3_vector *t_dest, su3_vector **ritz, Real *ritzVals,
			  quark_invert_control *qic, Real mass, int m, int nev, int deflation_grid, Real restart_tol,
			  imp_ferm_links_t *fn)
{

  char myname[] = "ks_eigcg_parity_gpu";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;//wrong.
#endif

//  if(qic->relresid != 0.){
//    printf("%s: GPU code does not yet support a Fermilab-type relative residual\n",myname);
//    terminate(1);
//  }
 
  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;

  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic->parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("eigcg: (right-hand side %d) source_norm = %e\n", curr_rhs_idx, (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY(i,qic->parity){
      memset(t_dest + i, 0, sizeof(su3_vector));
    } END_LOOP

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("EIGCG: time = %e (fn %s) iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

    return 0;
  }

  /* Initialize QUDA parameters */

  printf("Calling qudaEigCGInvert\n");
  fflush(stdout);

  initialize_quda();
 
  if(qic->parity == EVEN){
	  inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic->parity == ODD){
	  inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED) || defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif
  inv_args.mixed_precision = 0;//enforce full precision, for test only!!

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;
  const int ritz_precision = (inv_args.mixed_precision == 0) ? quda_precision : 1;//1 "means" single precision...

  double residual, relative_residual;
  int num_iters;
//WARNING: changed order!
  qudaEigCGInvert(PRECISION,
	          quda_precision, 
	          mass,
	          inv_args,
	          qic->resid,
	          qic->relresid,
	          fatlink, 
	          longlink,
                  u0,
	          t_src, 
	          t_dest,
                  ritz,//array of ritz vectors
                  ritzVals,
                  ritz_precision,
                  m, nev,
                  deflation_grid,
                  restart_tol,//e.g.: 5e+3*target_residual
                  curr_rhs_idx,//current rhs
                  last_rhs_flag,//is this the last rhs to solve?
	          &residual,
	          &relative_residual, 
	          &num_iters);
  

  qic->final_rsq    = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters  = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0; 

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("EIGCG: time = %e (fn_QUDA %s) iters = %d mflops = %e\n",
	   dtimec, prec_label[quda_precision-1], qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return num_iters;
}


int ks_eigcg_field( int curr_rhs, int tot_rhs_num, su3_vector *src, su3_vector *dest, su3_vector **ritzVectors, Real *ritzVals,
		      quark_invert_control *qic, Real mass, int m, int nev, int dgrid, Real restart_tol,  imp_ferm_links_t *fn)
{

  if(this_node==0){ printf("\nRunning RHS %d, parity = %d\n", curr_rhs, qic->parity);}
  
  int iters = 0;
  int parity = qic->parity;

  int last_rhs = curr_rhs < tot_rhs_num ? 0 : 1;

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_eigcg_parity_gpu(curr_rhs, last_rhs, src, dest, ritzVectors, ritzVals, qic,  mass, m, nev, dgrid, restart_tol, fn);
    report_status(qic);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_eigcg_parity_gpu(curr_rhs, last_rhs, src, dest, ritzVectors, ritzVals, qic,  mass, m, nev, dgrid, restart_tol, fn);
    report_status(qic);
  }

  qic->parity = parity;
  return iters;
}



int ks_eigcg_uml_field( int curr_rhs, int tot_rhs_num, su3_vector *src, su3_vector *dst, su3_vector **ritzVectors, Real *ritzVals,
                      quark_invert_control *qic, Real mass, int m, int nev, int dgrid, Real restart_tol,  imp_ferm_links_t *fn)

//int ks_eigcg_uml_field(su3_vector *src, su3_vector *dst,
//                         quark_invert_control *qic,
//                         Real mass, imp_ferm_links_t *fn )
{
    int cgn;
    register int i;
    register site *s;
    su3_vector *tmp = create_v_field();
    su3_vector *ttt = create_v_field();
    int even_iters;

    /* "Precondition" both even and odd sites */
    /* temp <- M_adj * src */

    ks_dirac_adj_op( src, tmp, mass, EVENANDODD, fn );

    /* dst_e <- (M_adj M)^-1 tmp_e  (even sites only) */
    qic->parity     = EVEN;
    cgn = ks_eigcg_field( curr_rhs, tot_rhs_num, src, dst, ritzVectors, ritzVals, qic,  mass, m, nev, dgrid, restart_tol, fn );
    even_iters = qic->final_iters;

    /* reconstruct odd site solution */
    /* dst_o <-  1/2m (Dslash_oe*dst_e + src_o) */
#ifdef FN
    dslash_fn_field( dst, ttt, ODD, fn );
#else
    //dslash_eo_field();
    //not applicable here: see dslash_ks_redefine.h
#endif
    FORODDSITES(i,s){
      sub_su3_vector( src+i, ttt+i, dst+i);
      scalar_mult_su3_vector( dst+i, 1.0/(2.0*mass), dst+i );
    }

    /* Polish off odd sites to correct for possible roundoff error */
    /* dst_o <- (M_adj M)^-1 temp_o  (odd sites only) */
    qic->parity = ODD;
    //!!!cgn += ks_congrad_field( tmp, dst, qic, mass, fn ); // or call init cg? think about this!

    //cgn = ks_eigcg_field( curr_rhs, tot_rhs, src, dest, ritzVectors, qic,  mass, m, nev, dgrid, restart_tol, fn );
    qic->final_iters += even_iters;

    //    check_invert_field( dst, src, mass, 1e-6, fn);
    destroy_v_field(tmp);
    destroy_v_field(ttt);
    
    return cgn;
}


