/********** update_ora.c ****************************************************/
/* MIMD version 6 */

/*
 Update lattice with microcanonical overrelaxed algorithm
*/

#include "generic_pg_includes.h"

int update()  {
int iters=0;

    /* check unitarity before doing anything */
    check_unitarity();
	
    /* do "steps" overrelaxed steps and stepsQ  qhb steps */
	relax(steps); 
	monte(stepsQ);

        /* reunitarize the gauge field */

        reunitarize();         

    if(steps > 0)return (iters/steps);
    else return(-99);
}
