/********** update.c ****************************************************/
/* MIMD version 7 */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

/*
 Update lattice.
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/

#include "symanzik_sl32_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>   /* For "finite" */
#endif

int update()  {
int step, iters=0;
#ifdef HMC_ALGORITHM
double startaction,endaction,change;
Real xrandom;
#endif

    /* refresh the momenta */
    ranmom();

    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){

#ifdef HMC_ALGORITHM
	/* find action */
	if(step==1){
	    startaction=d_action();
	    /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
	}
#endif

	/* update U's to middle of interval */
	if(step==1)update_u(0.5*epsilon);

	/* now update H by full time interval */
	update_h(epsilon);

	/* update U's by half time step to get to even time */
	if(step != steps)update_u(epsilon);
	else             update_u(0.5*epsilon);

	/* reunitarize the gauge field */
	reunitarize();

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
