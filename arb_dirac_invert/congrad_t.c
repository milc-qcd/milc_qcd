/******* congrad_t.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 4 */
/* Calls delta0.c for fermion matrix */

/* version of 21 Jan 99  */

/* This is an implementation of the conjugate gradient
algorithm a la version4/congrd_w.c */

/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector,
   "p", "mp", "tmp",  and "v" are working vectors for the
   conjugate gradient.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/


#include "arb_dirac_inv_includes.h"
#include <string.h>

int congrad_t( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    )
{
    int N_iter;
    register int i;
    register site *s;

  int iteration; /* counter for iterations */
  double rsq,source_norm,rsqmin,rsqstop;

    Real cp;
    Real test_tr;
    Real a, b, c, d;

  rsq = source_norm = 0.0;
  rsqmin=RsdCG*RsdCG;
  iteration=0;



    /* code if you want to start with psi=0... 
Kluge to work around Wilson code with point source*/
    if(start_flag == 0) {
        if(this_node==0)printf("psi_0=0\n");
        FORALLSITES(i,s) {
            clear_wvec( ((wilson_vector *)F_PT(s,dest)));
            copy_wvec(((wilson_vector *)F_PT(s,src)), 
			((wilson_vector *)F_PT(s,dest)));
        }
     }   


  delta0( dest, F_OFFSET(mp), PLUS );
  iteration++ ;	/* iteration counts number of multiplications
		   both by M and M_adjoint */
  total_iters++;




  FORALLSITES(i,s){
    sub_wilson_vector( ((wilson_vector *)F_PT(s,src)), &(s->mp), &(s->r) );
    source_norm += (double)magsq_wvec(((wilson_vector *)F_PT(s,src))  );
    rsq += (double)magsq_wvec( &(s->r) );
  }

  g_doublesum( &source_norm );
  g_doublesum( &rsq );

  if(this_node==0){
    printf("congrad: source_norm = %e\n",source_norm);
    fflush(stdout);
  } 

  rsqstop = rsqmin * source_norm;


    /*  p = M_dag*r  */
    delta0(F_OFFSET(r),F_OFFSET(p),MINUS);
    iteration++;
    total_iters++;



    /* cp = |p|^2 */
    cp=0.0;
    FORALLSITES(i,s) {
        cp += magsq_wvec( &(s->p) );
    }
    g_floatsum(&cp);
    /* test 
    if(this_node==0)    printf("trace of p--cp=%e\n",(double)cp);
     test */
    for( N_iter = 0; N_iter < MaxCG && rsq>rsqstop; 
        N_iter = N_iter + 1) {

        c=cp;
        /*   mp = M(u)*p */
        delta0(F_OFFSET(p),F_OFFSET(mp),PLUS);


        /* d = |mp|^2  */
        d=0.0;
        FORALLSITES(i,s) {
            d += magsq_wvec( &(s->mp) );
        }
        g_floatsum(&d);
        a = (c/d);
        /* test 
        if(this_node==0) printf(" top of congrad--c=%e, d=%e, a=%e\n",(double)c,
            (double)d,(double)a);
         test */

        /* psi = psi + a*p  */
        FORALLSITES(i,s) {
            scalar_mult_add_wvec(  ((wilson_vector *)F_PT(s,dest)),
                &(s->p), a, ((wilson_vector *)F_PT(s,dest))  );
        }

        /* r = r - a*mp */
        FORALLSITES(i,s) {
            scalar_mult_add_wvec( &(s->r),
                &(s->mp),-a, &(s->r) );
        }

        /*   mp = M(u)dag*r  */

        delta0(F_OFFSET(r),F_OFFSET(mp),MINUS);
    iteration++;
    total_iters++;


        /*   cp = |mp|^2  */
        cp=0.0;
        FORALLSITES(i,s) {
            cp += magsq_wvec( &(s->mp) );
        }
        g_floatsum(&cp);
        b = (cp/c);
        /*   p = mp + b*p  */
        FORALLSITES(i,s) {
            scalar_mult_add_wvec( &(s->mp),
                &(s->p), b, &(s->p) );
        }

        /* test 
        if(this_node==0) printf(" bottom of congrad--cp=%e, b=%e\n",(double)cp,
	    (double)b);
         test */

        rsq=0.0;
        FORALLSITES(i,s) {
	rsq += (double)magsq_wvec( &(s->r) );
        }
        g_doublesum(&rsq);

        if(this_node==0 && ((N_iter / 5)*5==N_iter) ){
	printf("iter %d residue %e\n",N_iter,
	    (double)(rsq));
	fflush(stdout);}
    }

    if( rsq > rsqstop ) {
        if(this_node==0)printf(" CONGRAD Not Converged\n");
    }

	*size_r=rsq;
    return(iteration);
}
