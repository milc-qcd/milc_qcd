/********** update_dense.c ****************************************************/
/* MIMD version 4 */

/*
 Update lattice with Metropolis algorithm, with "almost quenched"
 high density.
*/
#include "su3_dense_includes.h"

int update_dense()  {
   int Nhit;
   int step, iters=0;
   void relax_space(),monte_space(),monte_time(),heatbath_space();

    /* check unitarity before doing anything */
	check_unitarity();
	
    /* do "steps_over" overrelaxed steps and steps_update Metropolis steps */
	relax_space(steps_over); 
	monte_space(steps_update);
	/*heatbath_space(1);*/
	monte_time(steps_update);

        /* reunitarize the gauge field */

	if( phases_in != OFF ){
	  node0_printf("DUMMY: Reunitarizing with phases in!\n");
	  exit(0);
	}
        reunitarize();         

    if(steps_over > 0)return (iters/steps_over);
    else return(-99);
}
