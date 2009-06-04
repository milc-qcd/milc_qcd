/******* d_congrad2_cl.c - conjugate gradient for clover fermions ****/
/* MIMD version 7 */
/* Clover fermions */
/* For clover_dynamical/update.c.  Solves M_adjoint*M psi = chi */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
 
  without LU:
   M = A - kappa*( Dslash_eo + DSLASH_oe )
  with LU:
   M = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe
*/
#ifdef LU
#define FORMYSITES FOREVENSITESDOMAIN
#else
#define FORMYSITES FORALLSITESDOMAIN
#endif

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include "generic_clover_includes.h"

int congrad_cl(int niter,Real rsqmin,Real *final_rsq_ptr) 
{
register int i;
register site *s;
int iteration;	/* counter for iterations */
double source_norm;
double rsqstop;
Real a,b;
double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
void dslash_w_site_special();
msg_tag *tag[8],*tag2[8];
#ifdef LU
Real KAP = -kappa*kappa;
#else
Real KAP = -kappa;
#endif

double dtime;
dtime= -dclock();
	
	iteration=0;
start:
	/* mp <-  M_adjoint*M*psi
	   r,p <- chi - mp
	   rsq = |r|^2
	   source_norm = |chi|^2
	*/
	rsq = source_norm = 0.0;
#ifdef LU
	mult_this_ldu_site(gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVEN);
        dslash_w_site_special(F_OFFSET(psi), F_OFFSET(mp), PLUS, ODD, tag, 0);
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), ODD);
        dslash_w_site_special(F_OFFSET(tmp), F_OFFSET(mp), PLUS, EVEN, tag2, 0);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp));
        }
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), EVEN);
        dslash_w_site_special(F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD, tag, 1);
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), ODD);
        dslash_w_site_special(F_OFFSET(tmp), F_OFFSET(mp), MINUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
	    s->p = s->r;
	    source_norm += (double)magsq_wvec( &(s->chi) );
	    rsq += (double)magsq_wvec( &(s->r) );
        }
#else
	mult_this_ldu_site(gen_clov, F_OFFSET(psi), F_OFFSET(tmp), EVENANDODD);
	dslash_w_site_special(F_OFFSET(psi), F_OFFSET(mp), PLUS, EVENANDODD, tag, 0);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
	}
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), EVENANDODD);
	dslash_w_site_special(F_OFFSET(mp), F_OFFSET(mp), MINUS, EVENANDODD, tag, 1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
	    s->p = s->r;
	    source_norm += (double)magsq_wvec( &(s->chi) );
	    rsq += (double)magsq_wvec( &(s->r) );
	}
#endif
	g_doublesum( &source_norm );
	g_doublesum( &rsq );
        iteration++ ;	/* iteration counts number of multiplications
			   by M_adjoint*M */
	total_iters++;
/**if(this_node==0)printf("congrad2: source_norm = %e\n",source_norm);
if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
iteration,(double)rsq,(double)pkp,(double)a );**/

	rsqstop = rsqmin * source_norm;
	if( rsq <= rsqstop ){
	    *final_rsq_ptr= (Real)rsq;
	    for( i=XUP; i <= TUP; i++) {
		cleanup_gather(tag[i]);
		cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
		cleanup_gather(tag2[i]);
		cleanup_gather(tag2[OPP_DIR(i)]);
#endif
	    }
            cleanup_dslash_wtemps();
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
	mult_this_ldu_site(gen_clov, F_OFFSET(p), F_OFFSET(tmp), EVEN);
        dslash_w_site_special(F_OFFSET(p), F_OFFSET(mp), PLUS, ODD, tag, 1);
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), ODD);
        dslash_w_site_special(F_OFFSET(tmp), F_OFFSET(mp), PLUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
        }
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), EVEN);
        dslash_w_site_special(F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD, tag, 1);
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), ODD);
        dslash_w_site_special(F_OFFSET(tmp), F_OFFSET(mp), MINUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
            pkp += (double)wvec_rdot( &(s->p), &(s->mp) );
        }
#else
	mult_this_ldu_site(gen_clov, F_OFFSET(p), F_OFFSET(tmp), EVENANDODD);
	dslash_w_site_special(F_OFFSET(p), F_OFFSET(mp), PLUS, EVENANDODD, tag, 1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
	}
	mult_this_ldu_site(gen_clov, F_OFFSET(mp), F_OFFSET(tmp), EVENANDODD);
	dslash_w_site_special(F_OFFSET(mp), F_OFFSET(mp), MINUS, EVENANDODD, tag, 1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->tmp), &(s->mp), KAP, &(s->mp) );
            pkp += (double)wvec_rdot( &(s->p), &(s->mp) );
	}
#endif
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (Real)(rsq/pkp);
	rsq = 0.0;
	FORMYSITES(i,s){
            scalar_mult_add_wvec( &(s->psi), &(s->p), a, &(s->psi) );
            scalar_mult_add_wvec( &(s->r), &(s->mp), -a, &(s->r) );
	    rsq += (double)magsq_wvec( &(s->r) );
        }
	g_doublesum( &rsq );
/**if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
iteration,(double)rsq,(double)pkp,(double)a );**/

	if( rsq <= rsqstop ){
	    *final_rsq_ptr= (Real)rsq;
	    for( i=XUP; i <= TUP; i++) {
		cleanup_gather(tag[i]);
		cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
		cleanup_gather(tag2[i]);
		cleanup_gather(tag2[OPP_DIR(i)]);
#endif
	    }
dtime += dclock();
/** UMH if(this_node==0)printf("CONGRAD2: time = %e iters = %d mflops = %e\n",
dtime,iteration,(double)(2840.0*volume*iteration/(1.0e6*dtime*numnodes())) ); **/
             cleanup_dslash_wtemps();
	     return (iteration);
	}

	b = (Real)(rsq/oldrsq);
	FORMYSITES(i,s){
	    scalar_mult_add_wvec( &(s->r), &(s->p), b, &(s->p) );
	}

    } while( iteration%niter != 0);

    for( i=XUP; i <= TUP; i++) {
	cleanup_gather(tag[i]);
	cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
	cleanup_gather(tag2[i]);
	cleanup_gather(tag2[OPP_DIR(i)]);
#endif
    }
    if( iteration < 3*niter ) goto start;
    *final_rsq_ptr= (Real)rsq;
    if( rsq > rsqstop ){
	if(this_node==0)printf("No convergence in d_congrad2\n");
    }
    cleanup_dslash_wtemps();
    return(iteration);
}
