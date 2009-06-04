/********** update.c *** for clover fermions *************************/
/* MIMD version 7 */

/*
 Update lattice.
 Improved method for 1-2 flavors:
	update U by (epsilon/2)*(1-Nf/2)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/

#include "cl_dyn_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>    /* For "finite" */
#endif

void predict_next_psi(Real *oldtime,Real *newtime,Real *nexttime);

int update()  {
int step, iters=0;
Real final_rsq;
Real cg_time; /* simulation time for last two CG's */
#ifdef PHI_ALGORITHM
 Real old_cg_time,next_cg_time;
 double starttrlogA, endtrlogA;
#endif
Real CKU0 = kappa*clov_c/(u0*u0*u0);
#ifdef HMC_ALGORITHM
double startaction = 0, endaction, change;
Real xrandom;
#endif

/* printf("in update.c CKU0 = %f\n", CKU0); */
    /* refresh the momenta */
    ranmom();

    boundary_flip(MINUS);	/* turn on antiperiodic b.c. */
    /* do "steps" microcanonical steps"  */
    for(step=1; step <= steps; step++){

#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/

        if(step==1){
          compute_clov(gen_clov, CKU0);
#ifdef LU
	  compute_clovinv(gen_clov, ODD);
	  starttrlogA = gen_clov->trlogA;
#else
	    starttrlogA = (double)0.0;
#endif /*LU*/
	    grsource_w();
	    old_cg_time = cg_time = -1.0e6;
 	}

#ifdef HMC_ALGORITHM
        /* find action */
        if(step==1){
          /* do conjugate gradient to get (Madj M)inverse * chi */
	  iters += congrad_cl(niter,rsqmin,&final_rsq);
	  cg_time = 0.0;
        /**checkmul();**/
     	   startaction=d_action();
     	   startaction -= (double)2.0 * starttrlogA;
	   /* printf("startaction= %g\n",startaction); */
           /* copy link field to old_link */
	   gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
        }
#endif /*hmc*/
        if(step==1){
	  free_this_clov(gen_clov);
	}
	/* update U's to middle of interval */
     	update_u(0.5*epsilon);

        /* save conjugate gradient solution, predict next one */
        next_cg_time = ((Real)step-0.5)*epsilon;
        predict_next_psi(&old_cg_time,&cg_time,&next_cg_time);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
       	update_u(epsilon*(0.5-nflavors/4.0));

        /* generate a pseudofermion configuration */
	compute_clov(gen_clov, CKU0);
#ifdef LU
	compute_clovinv(gen_clov, ODD);
#endif /*LU*/
     	grsource_w();
	free_this_clov(gen_clov);

	/* update U's to middle of interval */
     	update_u(epsilon*nflavors/4.0);
#endif /* phi & R */

        /* do conjugate gradient to get (Madj M)inverse * chi */
	compute_clov(gen_clov, CKU0);
#ifdef LU
	compute_clovinv(gen_clov, ODD);
#endif /*LU*/
	  iters += congrad_cl(niter,rsqmin,&final_rsq);
        cg_time = ((Real)step - 0.5)*epsilon;
/**checkmul();**/

	/* now update H by full time interval */
    	update_h(epsilon);
	free_this_clov(gen_clov);

    	/* update U's by half time step to get to even time */
    	update_u(epsilon*0.5);

        /* reunitarize the gauge field */
	boundary_flip(PLUS);
        reunitarize();
	boundary_flip(MINUS);

    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * chi */
    next_cg_time = steps*epsilon;
    predict_next_psi(&old_cg_time,&cg_time,&next_cg_time);
    compute_clov(gen_clov, CKU0);
#ifdef LU
    compute_clovinv(gen_clov, ODD);
    endtrlogA = gen_clov->trlogA;
#else
    endtrlogA = (double)0.0;
#endif /*LU*/
    iters += congrad_cl(niter,rsqmin,&final_rsq);
    free_this_clov(gen_clov);
    cg_time = steps*epsilon;
    endaction=d_action();
    endaction -= (double)2.0 * endtrlogA;
    /* printf("endaction= %g\n",endaction); */
    change = endaction-startaction;
    /* Reject configurations giving overflow */
#ifndef HAVE_IEEEFP_H
    if(fabs((double)change)>1e20){
#else
    if(!finite((double)change)){
#endif
	if(this_node==0)printf(
	    "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
	change = 1.0e20;
    }

    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( -change ) < (double)xrandom ){
	if(steps > 0)
	    gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
	if(this_node==0)printf("REJECT: delta S = %e\n", change);
    }
    else {
	if(this_node==0)printf("ACCEPT: delta S = %e\n", change);
    }
#endif /*HMC*/

    boundary_flip(PLUS);
    if(steps > 0)return (iters/steps);
    else return(-99);
}

#ifdef PHI_ALGORITHM
#ifdef LU
#define FORMYSITES FOREVENSITES
#else
#define FORMYSITES FORALLSITES
#endif
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
void predict_next_psi(Real *oldtime,Real *newtime,Real *nexttime)
{
register int i;
register site *s;
register Real x;
wilson_vector tvec;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
	FORMYSITES(i,s){
	    s->old_psi = s->psi;
	}
    }
    else  {
	FORMYSITES(i,s){
            sub_wilson_vector( &(s->psi), &(s->old_psi), &tvec);
	    s->old_psi = s->psi;
            scalar_mult_add_wvec( &(s->psi), &tvec,x, &(s->psi) );
	}
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}
#endif

