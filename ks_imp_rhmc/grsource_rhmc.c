/************************ grsource_rhmc.c *****************************/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*M)^(nf/8)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

// residues, roots and order define rational function approximation for
 // x^(nf/8)
void grsource_imp_rhmc( field_offset dest, params_ratfunc *rf,
			int parity, su3_vector **multi_x, su3_vector *sumvec,
			Real my_rsqmin, int my_niter, int my_prec,
			imp_ferm_links_t *fn, int naik_term_epsilon_index,
			Real naik_term_epsilon)
{
  register int i,j;
  register site *s;
  Real final_rsq;
  int order = rf->order;
  Real *residues = rf->res;
  Real *roots = rf->pole;
  /*TEMP*/ double sum;
  double dtimec = -dclock();
  
  sum=0.0;
  FORSOMEPARITY_OMP(i,s,parity,private(j) reduction(+:sum)){
    for(j=0;j<3;j++){
#ifdef SITERAND
      s->g_rand.c[j] = complex_gaussian_rand_no(&(s->site_prn));
#else
      s->g_rand.c[j] = complex_gaussian_rand_no(&node_prn);
#endif
    }
    /*TEMP*/ sum += (double)magsq_su3vec( &(s->g_rand) );
  } END_LOOP_OMP
  /*TEMP*/g_doublesum( &sum);  node0_printf("GRSOURCE: sum = %.10e\n",sum);
  dtimec += dclock();
  ks_ratinv( F_OFFSET(g_rand), multi_x, roots, residues, order, my_niter,
	     my_rsqmin, my_prec, parity, &final_rsq, fn, 
	     naik_term_epsilon_index, naik_term_epsilon );
  dtimec -= dclock();
  ks_rateval( sumvec, F_OFFSET(g_rand), multi_x, residues, order, parity );
  //  dtimec -= dclock();
  FORSOMEPARITY_OMP(i,s,parity,default(shared) ){ *(su3_vector *)F_PT(s,dest) = sumvec[i]; } END_LOOP_OMP
  dtimec += dclock();
  node0_printf("GRSOURCETIME: time = %e\n",dtimec);
}/* grsource_rhmc */
