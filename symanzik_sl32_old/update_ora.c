/********** update_ora.c ****************************************************/
/* MIMD version 7 */

/*
 Update lattice with microcanonical overrelaxed algorithm
*/

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#include "symanzik_sl32_includes.h"

int update()  {
int iters=0;

    /* do "steps" overrelaxed steps and stepsQ  qhb steps */
	relax(steps); 
	monte(stepsQ);

	/* reunitarize the gauge field */
	reunitarize();         

    if(steps > 0)return (iters/steps);
    else return(-99);
}
