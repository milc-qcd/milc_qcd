/********** update.c ****************************************************/
/* MIMD version 6 */

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

#include "wi_dyn_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>      /* For "finite" */
#endif

int update()  {
int step, iters=0;
Real final_rsq;
void update_u( Real eps ), update_h( Real eps);
void predict_next_psi(Real *oldtime,Real *newtime,Real *nexttime);
Real cg_time,old_cg_time,next_cg_time; /* simulation time for last two CG's */
#ifdef HMC_ALGORITHM
double startaction,endaction,change,d_action();
Real xrandom;
#endif

    /* refresh the momenta */
    ranmom();

    boundary_flip(MINUS);	/* turn on antiperiodic b.c. */
    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){

#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/
        if(step==1){grsource_w(); old_cg_time = cg_time = -1.0e6;}

#ifdef HMC_ALGORITHM
        /* find action */
        /* do conjugate gradient to get (Madj M)inverse * chi */
        if(step==1){
            /* do conjugate gradient to get (Madj M)inverse * chi */
	    iters += congrad_w( niter,rsqmin,&final_rsq);
	    cg_time = 0.0;
/**checkmul();**/
     	    startaction=d_action();
            /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
        }
#endif

	/* update U's to middle of interval */
     	update_u(0.5*epsilon);

        /* save conjugate gradient solution, predict next one */
        next_cg_time = ((Real)step-0.5)*epsilon;
        predict_next_psi(&old_cg_time,&cg_time,&next_cg_time);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
       	update_u(epsilon*(0.5-nflavors/4.0));

        /* generate a pseudofermion configuration */
     	grsource_w();

	/* update U's to middle of interval */
     	update_u(epsilon*nflavors/4.0);
#endif

        /* do conjugate gradient to get (Madj M)inverse * chi */
	iters += congrad_w( niter,rsqmin,&final_rsq);
        cg_time = ((Real)step - 0.5)*epsilon;
/**checkmul();**/

	/* now update H by full time interval */
    	update_h(epsilon);

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
    iters += congrad_w( niter,rsqmin,&final_rsq);
    cg_time = steps*epsilon;
    endaction=d_action();
    change = endaction-startaction;
#ifndef HAVE_IEEEFP_H
    /* Floating overflow gives +- NaN for change - need to reject if -NaN */
    if(change < -1.0e20 ){
	if(this_node==0)printf(
	    "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
	change = 1.0e20;
    }
#else
    /* Carleton's change for detecting NaN */
    /* Reject configurations giving overflows */
    if(!finite((double)change)){
	if(this_node==0)printf(
	    "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
	change = 1.0e20;
    }
#endif

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
#endif

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
