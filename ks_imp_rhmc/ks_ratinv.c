/******* ks_ratinv.c - multi-mass CG for rational function approximations
 for SU3/fermions 
 Also includes a routine to evaluate the rational function after the inversion
 is done
****/
/* MIMD version 7 */

/* Multi-mass CG inverter for staggered fermions */

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. 

   This version is based on d_congrad5_fn.c and d_congrad5_eo.c 

   approximation form:  a0 + a1/(x-b1) + ... aN/(x-bN)  t
   N = "order"
   a[i] = residues
   b[i] = roots
   psim[i-1] = 1/(M^dagger M + b[i] )

   For "fat link actions", ie when FN is defined, this version
   assumes connection to nearest neighbor points is stored in fatlink.
   For actions with a Naik term, it assumes the connection to third
   nearest neighbors is in longlink.

*/


#include "ks_imp_includes.h"	/* definitions files and prototypes */

int ks_ratinv(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *roots,	/* the roots */
    int order,		/* order of rational function approx */
    int my_niter,	/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* desired intermediate precicion */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    imp_ferm_links_t *fn_const, /* Fermion links */
    int naik_term_epsilon_index, /* Index of naik term common to this set */
    Real naik_term_epsilon /* Epsilon common to this set */
    )
{
    // Just a multimass inversion.  start at roots[1] because first term
    // in rational function expansion is just the constant term
    // It's "order" instead of "order-1" because order is order of expansion, arrays
    // have order+1 elements 
  quark_invert_control *qic;
  su3_vector *in;
  ks_param *ksp;
  int k;
  imp_ferm_links_t **fn;
  char myname[] = "update_h_rhmc";
  int iters;

  in = create_v_field_from_site_member(src);

  /* Set up inversion control structure */
  /* For molecular dynamics they are identical */
  qic = (quark_invert_control *)malloc(order*sizeof(quark_invert_control));
  if(qic == NULL){
    printf("ks_ratinv: No room for qic\n");
    terminate(1);
  }

  for(k = 0; k < order; k++){
    qic[k].prec = prec;
    qic[k].min = 0;
    qic[k].max = my_niter;
    qic[k].nrestart = nrestart;
    qic[k].parity = parity;
    qic[k].start_flag = 0;
    qic[k].nsrc = 1;
    qic[k].resid = sqrt(rsqmin);
    qic[k].relresid = 0;
  }

  /* Load ks parameters for inverters */
  ksp = (ks_param *)malloc(order*sizeof(ks_param));
  if(ksp == NULL){
    printf("%s(%d): No room\n", myname, this_node);
    terminate(1);
  }
  
  for(k = 0; k < order; k++){
    ksp[k].offset = roots[k+1];
#if FERM_ACTION == HISQ
    ksp[k].naik_term_epsilon = naik_term_epsilon;
    ksp[k].naik_term_epsilon_index = naik_term_epsilon_index;
#endif
  }

  /* Set fn links for the inversion (all the same here) */
  fn = (imp_ferm_links_t **)malloc(order*sizeof(imp_ferm_links_t *));
  if(fn == NULL){
    printf("%s(%d): No room\n", myname, this_node);
    terminate(1);
  }
  
  for(k = 0; k < order; k++)
    fn[k] = fn_const;

  iters = ks_multicg_field( in, psim, ksp, order, qic, fn );

  free(fn);
  free(ksp);

  destroy_v_field(in);
  *final_rsq_ptr = qic[0].final_rsq;
  free(qic);
  return iters;
}

/* evaluate the rational function approximation after all the
  1 / ( M^dagger M + root) 's have been computed */
int ks_rateval(	
    su3_vector *dest,	/* answer vector */
    field_offset src,	/* source vector (for a_0 term) */
    su3_vector **psim,	/* solution vectors  from multiroot CG */
    Real *residues,	/* the residues */
    int order,		/* order of approximation */
    int parity		/* parity to be worked on */
    )
{
   register int i,j; register site *s;
   FORSOMEPARITY(i,s,parity){
      scalar_mult_su3_vector( (su3_vector *)F_PT(s,src), residues[0], &(dest[i]) );
      for(j=1;j<=order;j++){
        scalar_mult_add_su3_vector( &(dest[i]), &(psim[j-1][i]), residues[j], &(dest[i]) );
      }
   }
   return 0;
}
