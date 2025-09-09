/*************** d_action_rhmc.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.  */

#include <ks_imp_includes.h>	/* definitions files and prototypes */
#include "../include/fermion_links.h"
#include "../include/openmp_defs.h"
#include "../include/generic_quda.h"

static Real ahmat_mag_sq(anti_hermitmat *pt);

/*DEBUG*/
double old_g, old_h, old_f, old_a;
/*ENDDEBUG*/

double d_action_rhmc( su3_vector **multi_x, su3_vector *sumvec){
  double ssplaq,stplaq,g_action,h_action,f_action;
  double dtimec = -dclock();
  
  plaquette_action(&ssplaq,&stplaq);
  ssplaq *= -1.0; stplaq *= -1.0;
  g_action = -beta*volume*(ssplaq+stplaq);
  node0_printf("PLAQUETTE ACTION: %e\n",g_action);
  
  g_action = (beta/3.0)*imp_gauge_action_ks();
  h_action = hmom_action();
  f_action = fermion_action(multi_x,sumvec);
  
  node0_printf("ACTION: g,h,f = %.14e  %.14e  %.14e  %.14e\n",
	       g_action, h_action, f_action, g_action+h_action+f_action );
  
  /*DEBUG*/
  node0_printf("DG = %e, DH = %e, DF = %e, D = %e\n",
	       g_action-old_g, h_action-old_h, f_action-old_f,
	       g_action+h_action+f_action-old_a);
  old_g=g_action; old_h=h_action; old_f=f_action;
  old_a=g_action+h_action+f_action;
  /*ENDDEBUG*/
  
  dtimec += dclock();
  node0_printf("ACTIONTIME: time = %e\n",dtimec);
  return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
double fermion_action( su3_vector **multi_x, su3_vector *sumvec) {
  register int i;
  register site *s;
  Real final_rsq;
  double sum;
  int iphi, inaik, jphi, n; 
  imp_ferm_links_t *fn;
  sum=0.0;
  iphi=0;
#if FERM_ACTION == HISQ
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( inaik=0; inaik<n; inaik++ ) {
    restore_fermion_links_from_site(fn_links, prec_fa[0]);
    fn = get_fm_links(fn_links, inaik);
    for( jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
      ks_ratinv( F_OFFSET(phi[iphi]), multi_x,
                 rparam[iphi].FA.pole, rparam[iphi].FA.res,
		 rparam[iphi].FA.order, niter_fa[iphi], rsqmin_fa[iphi], 
		 prec_fa[iphi], EVEN, &final_rsq, fn, 
		 inaik, rparam[iphi].naik_term_epsilon );
      ks_rateval( sumvec, F_OFFSET(phi[iphi]), multi_x, 
  		rparam[iphi].FA.res, rparam[iphi].FA.order, EVEN );
      FOREVENFIELDSITES_OMP(i,default(shared) reduction(+:sum)){ /* phi is defined on even sites only */
        sum += magsq_su3vec( &(sumvec[i]) );
      } END_LOOP_OMP
      iphi++;
    }
    destroy_fn_links(fn);
  }
  
  g_doublesum( &sum );
  return(sum);
}

/* plaquette action */
void plaquette_action(double *ss_plaq, double *st_plaq) {
#if !defined(HAVE_QUDA) || defined(SCHROED_FUN)
  d_plaquette(ss_plaq, st_plaq);
#else
  d_plaquette_gpu(ss_plaq, st_plaq);
  /* Apply a negative sign to "cancel" the negative sign above
     The CPU codepath assumes the phases are still on the field,
     which adds another minus sign to the plaquette */
  *ss_plaq = -*ss_plaq;
  *st_plaq = -*st_plaq;
#endif
}

/* gauge momentum contribution to the action */
double hmom_action(void) {
#ifndef HAVE_QUDA
  register int i,dir;
  register site *s;
  double sum;

  sum=0.0;
  FORALLSITES_OMP(i,s,private(dir) reduction(+:sum)) {
    for(dir=XUP;dir<=TUP;dir++){
      sum += (double)ahmat_mag_sq( &(s->mom[dir]) ) - 4.0;
      /* subtract 1/2 per d.o.f. to help numerical acc. in sum */
    }
  } END_LOOP_OMP
  g_doublesum( &sum );

  return(sum);
#else
  QudaMILCSiteArg_t arg = newQudaMILCSiteArg();
  return qudaMomAction(MILC_PRECISION, &arg);
#endif
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt){
  register Real x,sum;
  x = pt->m00im; sum  = 0.5*x*x;
  x = pt->m11im; sum += 0.5*x*x;
  x = pt->m22im; sum += 0.5*x*x;
  x = pt->m01.real; sum += x*x;
  x = pt->m01.imag; sum += x*x;
  x = pt->m02.real; sum += x*x;
  x = pt->m02.imag; sum += x*x;
  x = pt->m12.real; sum += x*x;
  x = pt->m12.imag; sum += x*x;
  return(sum);
}

