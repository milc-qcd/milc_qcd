/************************ grsource_rhmc.c *****************************/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"

/* construct a gaussian random vector, g_rand, and phi=M(dagger)*M)^(nf/8)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

// residues, roots and order define rational function approximation for
 // x^(nf/8)

/* Here the dest is in the site structure */
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
  } END_LOOP_OMP;
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
} /* grsource_imp_rhmc */

/* This version takes dest to be a field */
void grsource_imp_rhmc_field( su3_vector *dest, params_ratfunc *ratf,
			      int parity, su3_vector **multi_x, su3_vector *sumvec,
			      Real my_rsqmin, int my_niter, int my_prec,
			      imp_ferm_links_t *fn, int naik_term_epsilon_index,
			      Real naik_term_epsilon)
{
  register int i,j;
  register site *s;
  Real final_rsq;
  int order = ratf->order;
  Real *residues = ratf->res;
  Real *roots = ratf->pole;
  /*TEMP*/ double sum;
  double dtimec = -dclock();
  
  su3_vector *g_rand = create_v_field();
  sum=0.0;
  FORSOMEPARITY_OMP(i,s,parity,private(j) reduction(+:sum)){
    for(j=0;j<3;j++){
#ifdef SITERAND
      g_rand[i].c[j] = complex_gaussian_rand_no(&(s->site_prn));
#else
      g_rand[i].c[j] = complex_gaussian_rand_no(&node_prn);
#endif
    }
    /*TEMP*/ sum += (double)magsq_su3vec( &(g_rand[i]) );
  } END_LOOP_OMP
  /*TEMP*/g_doublesum( &sum);  node0_printf("GRSOURCE: sum = %.10e\n",sum);
  dtimec += dclock();
  ks_ratinv_field( g_rand, multi_x, roots, residues, order, my_niter,
		   my_rsqmin, my_prec, parity, &final_rsq, fn, 
		   naik_term_epsilon_index, naik_term_epsilon );
  dtimec -= dclock();
  ks_rateval_field( sumvec, g_rand, multi_x, residues, order, parity );
  FORSOMEPARITY_OMP(i,s,parity,default(shared) ){
    su3vec_copy(sumvec+i, dest+i);
  }
  END_LOOP_OMP;
  dtimec += dclock();
  node0_printf("GRSOURCETIME: time = %e\n",dtimec);

  destroy_v_field(g_rand);

} /* grsource_imp_rhmc_field */
