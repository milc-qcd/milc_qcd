/********** update_hgh.c *************************************************/
/* MIMD version 7 */
/* THIS CODE NEEDS UPGRADING NOW */
/*
 Update lattice.
 Momentum, then gauge, then momentum, ....
 but with the Sexton-Weingarten split of fermionic and pure gauge
	compute X
	update H to epsilon/2 with fermionic force only
	Do n_sxw times{
	    update H to epsilon/(2*n_sxw) with pure gauge force
	    update U to epsilon/n_sxw
	    update H to epsilon/(2*n_sxw) with pure gauge force
	}
	compute X
	update H full step with fermionic force only
	Do n_sxw times{
	    etc.

	compute X
	update H, to epsilon/2 with fermionic force only

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "schroed_ks_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>    /* For "finite" */
#endif

static void predict_next_xxx(Real *oldtime, Real *newtime, Real *nexttime);

int update()  {
int step, i_sxw, iters=0;
Real final_rsq, eps_sxw;
Real cg_time,old_cg_time,next_cg_time;	/* simulation time for last two CG's */
#ifdef HMC_ALGORITHM
double startaction,endaction,change;
Real xrandom;
#endif

    /* refresh the momenta */
    ranmom();

#ifdef SEXT_WEIN
    eps_sxw = epsilon/n_sxw;
#endif

    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){

#ifdef PHI_ALGORITHM
	/* generate a pseudofermion configuration only at start*/
	if(step==1){
	    grsource(EVEN);
	    old_cg_time = cg_time = -1.0e6;

	    /* do conjugate gradient to get (Madj M)inverse * phi */
	    /* NOTE: NEED TO UPGRADE TO ASQTAD.  BUILD ks_act_paths, ETC. */
	    load_ferm_links(&fn_links);
	    iters += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				niter, nrestart, rsqmin, PRECISION, EVEN, 
				&final_rsq, &fn_links);
	    cg_time = 0.0;
	}
#ifdef HMC_ALGORITHM
	/* find action */
	if(step==1){
	    startaction=d_action();
	    /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
	}
#endif

	if(step==1){
	    /* update H to middle of first interval */
#ifdef SEXT_WEIN
	    /* with fermion_force only! */
	    /* First compute M*xxx in temporary vector ttt */
	    /* The diagonal term in M doesn't matter */
	    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD );
	    fermion_force(0.5*epsilon);
#else
	    update_h(0.5*epsilon);
#endif
	}


#else /* "R" algorithm */
	/*Don't know how to do this yet */
	if(this_node==0)printf("OOPS: R algorithm not implemented\n");terminate();
#endif

#ifdef SEXT_WEIN
	/* update H by eps_sxw/2 with gauge force only */
	gauge_force(0.5*eps_sxw);
	for(i_sxw=1; i_sxw<=n_sxw; i_sxw++){
	    /* update U by eps_sxw */
	    update_u(eps_sxw);
	    if( i_sxw < n_sxw ){
		/* update H by eps_sxw with gauge force only */
		gauge_force(eps_sxw);
	    }
	    else{
		/* update H by eps_sxw/2 with gauge force only */
		gauge_force(0.5*eps_sxw);
	    }
	}
#else
	/* now update U by full time interval */
	update_u(epsilon);
#endif

	/* reunitarize the gauge field */
	rephase_sf( OFF );
	reunitarize();
	rephase_sf( ON );

	/* do conjugate gradient to get (Madj M)inverse * phi */
	next_cg_time = step*epsilon;
	predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
	load_ferm_links(&fn_links);
	iters += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, nrestart, rsqmin, PRECISION, EVEN, 
			    &final_rsq, &fn_links);
	cg_time = step*epsilon;

	if( step < steps ){
	    /* update H by full time step */
#ifdef SEXT_WEIN
	    /* with fermion_force only! */
	    /* First compute M*xxx in temporary vector ttt */
	    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD );
	    fermion_force(epsilon);
#else
	    update_h(epsilon);
#endif
	}
	else{
	    /* update H by half time step to get to even time */
#ifdef SEXT_WEIN
	    /* with fermion_force only! */
	    /* First compute M*xxx in temporary vector ttt */
	    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD );
	    fermion_force(0.5*epsilon);
#else
	    update_h(0.5*epsilon);
#endif
	}

    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
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

    if(steps > 0)return (iters/steps);
    else return(-99);
}

#ifdef PHI_ALGORITHM
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
static void predict_next_xxx(Real *oldtime, Real *newtime, Real *nexttime) {
register int i;
register site *s;
register Real x;
su3_vector tvec;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
        FOREVENSITES(i,s) if(s->t > 0){
	    s->old_xxx = s->xxx;
        }
    }
    else  {
        FOREVENSITES(i,s) if(s->t > 0){
            sub_su3_vector( &(s->xxx), &(s->old_xxx), &tvec);
	    s->old_xxx = s->xxx;
            scalar_mult_add_su3_vector( &(s->xxx), &tvec,x, &(s->xxx) );
        }
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}
#endif
