/******* congrad2.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 6 */
/* Wilson fermions */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
 
  without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
  with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
*/
#ifdef LU
#define FORMYSITES FOREVENSITES
#else
#define FORMYSITES FORALLSITES
#endif

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "phi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include "wi_dyn_includes.h"

int congrad(int niter,Real rsqmin,Real *final_rsq_ptr) 
{
register int i;
register site *s;
int iteration;	/* counter for iterations */
int iteration_mod_niter;
Real source_norm;
Real rsqstop;
Real a,b,rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
complex cc;
register int first,count;

/**double dtime,dclock();
dtime = -dclock();**/
	
	iteration=0;
start:
	/* mp <-  M_adjoint*M*psi
	   r,p <- chi - mp
	   rsq = |r|^2
	   source_norm = |chi|^2
	*/
	rsq = source_norm = 0.0;
#ifdef LU
        dslash(F_OFFSET(psi),F_OFFSET(psi) ,PLUS,ODD );
        dslash(F_OFFSET(psi),F_OFFSET(ttt),PLUS,EVEN);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->psi),&(s->ttt), -kappa*kappa,&(s->ttt) );
        }
        dslash(F_OFFSET(ttt),F_OFFSET(ttt),MINUS,ODD );
        dslash(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVEN);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa*kappa, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
/**scalar_mult_add_wvec( &(s->chi), &(s->mp), -1.0, &(s->r) );**/
	    s->p = s->r;
	    source_norm += magsq_wvec( &(s->chi) );
	    rsq += magsq_wvec( &(s->r) );
        }
#else
	dslash(F_OFFSET(psi),F_OFFSET(ttt),PLUS,EVENANDODD);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->psi), &(s->ttt), -kappa, &(s->ttt) );
	}
	dslash(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVENANDODD);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
/**scalar_mult_add_wvec( &(s->chi), &(s->mp), -1.0, &(s->r) );**/
	    s->p = s->r;
	    source_norm += magsq_wvec( &(s->chi) );
	    rsq += magsq_wvec( &(s->r) );
	}
#endif
	g_floatsum( &source_norm );
	g_floatsum( &rsq );
        iteration++ ;	/* iteration counts number of multiplications
			   by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	if( rsq <= rsqstop ){
	    *final_rsq_ptr=rsq;
	     return (iteration);
	}

    /* main loop - do until convergence or time to restart */
	/* 
	   oldrsq <- rsq
	   mp <- M_adjoint*M*p
	   pkp <- p.M_adjoint*M.p
	   a <- rsq/pkp
	   psi <- psi + a*p
	   r <- r - a*mp
	   rsq <- |r|^2
	   b <- rsq/oldrsq
	   p <- r + b*p
	*/
    do{
	oldrsq = rsq;
	pkp = 0.0;
#ifdef LU
        dslash(F_OFFSET(p),F_OFFSET(p) ,PLUS,ODD );
        dslash(F_OFFSET(p),F_OFFSET(ttt),PLUS,EVEN);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->p), &(s->ttt), -kappa*kappa, &(s->ttt) );
        }
        dslash(F_OFFSET(ttt),F_OFFSET(ttt),MINUS,ODD );
        dslash(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVEN);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa*kappa, &(s->mp) );
            pkp += wvec_rdot( &(s->p), &(s->mp) );
        }
#else
	dslash(F_OFFSET(p),F_OFFSET(ttt),PLUS,EVENANDODD);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->p), &(s->ttt), -kappa, &(s->ttt) );
	}
	dslash(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVENANDODD);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa, &(s->mp) );
            pkp += wvec_rdot( &(s->p), &(s->mp) );
	}
#endif
	g_floatsum( &pkp );
	iteration++;
	total_iters++;

	a = rsq/pkp;
	rsq = 0.0;
	FORMYSITES(i,s){
            scalar_mult_add_wvec( &(s->psi), &(s->p), a, &(s->psi) );
            scalar_mult_add_wvec( &(s->r), &(s->mp), -a, &(s->r) );
	    rsq += magsq_wvec( &(s->r) );
        }
	g_floatsum( &rsq );
/**if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
iteration,(double)rsq,(double)pkp,(double)a );**/

	if( rsq <= rsqstop ){
	    *final_rsq_ptr=rsq;
/**dtime += dclock();
if(this_node==0)printf("CONGRAD2: time = %e iters = %d mflops = %e\n",
dtime,iteration,(double)(2840.0*volume*iteration/(1.0e6*dtime*numnodes())) );**/
	     return (iteration);
	}

	b = rsq/oldrsq;
	FORMYSITES(i,s){
	    scalar_mult_add_wvec( &(s->r), &(s->p), b, &(s->p) );
	}

    } while( iteration%niter != 0);

    if( iteration < 2*niter ) goto start;
    *final_rsq_ptr=rsq;
    return(iteration);
}

