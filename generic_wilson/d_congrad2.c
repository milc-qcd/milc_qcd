/******* d_congrad2.c - conjugate gradient for SU3/fermions ****/
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

#include "generic_wilson_includes.h"

/**#define CGTIME**/

int congrad(int niter,Real rsqmin,Real *final_rsq_ptr) 
{
register int i;
register site *s;
int iteration;	/* counter for iterations */
double source_norm;
double rsqstop;
Real a,b;
double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
msg_tag *tag[8],*tag2[8];

double dtime;
#ifdef CGTIME
dtime= -dclock();
#endif
	
	iteration=0;
start:
	/* mp <-  M_adjoint*M*psi
	   r,p <- chi - mp
	   rsq = |r|^2
	   source_norm = |chi|^2
	*/
	rsq = source_norm = 0.0;
#ifdef LU
        dslash_special(F_OFFSET(psi),F_OFFSET(psi),PLUS,ODD, tag, 0 );
        dslash_special(F_OFFSET(psi),F_OFFSET(ttt),PLUS,EVEN,tag2,0 );
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->psi),&(s->ttt), -kappa*kappa,&(s->ttt) );
        }
        dslash_special(F_OFFSET(ttt),F_OFFSET(ttt),MINUS,ODD,tag,1 );
        dslash_special(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVEN,tag2,1 );
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa*kappa, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
/**scalar_mult_add_wvec( &(s->chi), &(s->mp), (Real)(-1.0), &(s->r) );**/
	    s->p = s->r;
	    source_norm += (double)magsq_wvec( &(s->chi) );
	    rsq += (double)magsq_wvec( &(s->r) );
        }
#else
	dslash_special(F_OFFSET(psi),F_OFFSET(ttt),PLUS,EVENANDODD,tag,0);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->psi), &(s->ttt), -kappa, &(s->ttt) );
	}
	dslash_special(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVENANDODD,tag,1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa, &(s->mp) );
	    sub_wilson_vector( &(s->chi), &(s->mp), &(s->r) );
/**scalar_mult_add_wvec( &(s->chi), &(s->mp), (Real)(-1.0), &(s->r) );**/
	    s->p = s->r;
	    source_norm += (double)magsq_wvec( &(s->chi) );
	    rsq += (double)magsq_wvec( &(s->r) );
	}
#endif
	g_doublesum( &source_norm );
/**{Real xxx; xxx=source_norm; g_floatsum( &xxx ); source_norm=xxx;}**/
	g_doublesum( &rsq );
/**{Real xxx; xxx=rsq; g_floatsum( &xxx ); rsq=xxx;}**/
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
        dslash_special(F_OFFSET(p),F_OFFSET(p) ,PLUS,ODD,tag,1 );
        dslash_special(F_OFFSET(p),F_OFFSET(ttt),PLUS,EVEN,tag2,1);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->p), &(s->ttt), -kappa*kappa, &(s->ttt) );
        }
        dslash_special(F_OFFSET(ttt),F_OFFSET(ttt),MINUS,ODD,tag,1 );
        dslash_special(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVEN,tag2,1);
        FOREVENSITES(i,s){
            scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa*kappa, &(s->mp) );
            pkp += (double)wvec_rdot( &(s->p), &(s->mp) );
        }
#else
	dslash_special(F_OFFSET(p),F_OFFSET(ttt),PLUS,EVENANDODD,tag,1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->p), &(s->ttt), -kappa, &(s->ttt) );
	}
	dslash_special(F_OFFSET(ttt),F_OFFSET(mp),MINUS,EVENANDODD,tag,1);
	FORALLSITES(i,s){
	    scalar_mult_add_wvec( &(s->ttt), &(s->mp), -kappa, &(s->mp) );
            pkp += (double)wvec_rdot( &(s->p), &(s->mp) );
	}
#endif
	g_doublesum( &pkp );
/**{Real xxx; xxx=pkp; g_floatsum( &xxx ); pkp=xxx;}**/
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
/**{Real xxx; xxx=rsq; g_floatsum( &xxx ); rsq=xxx;}**/
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
#ifdef CGTIME
dtime += dclock();
if(this_node==0)
  printf("CONGRAD2: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	 dtime,rsq,iteration,
	 (double)5616*iteration*even_sites_on_node/(dtime*1e6));
#endif
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
    return(iteration);
}

