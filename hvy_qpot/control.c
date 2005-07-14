/*********************  control.c ***********************/
/* MIMD version 6 */
/* Main procedure for pure gauge SU3 */

/* original code by UMH */
/* 2/19/98 Version 5 port CD */

/* This version does some required number of smearing iterations,
   then transforms to axial gauge and computes simple, i.e. one-
   plaquette, glueball operators and (time-like) Wilson loops
   for the computation of the heavy quark potential. */

#define CONTROL
#include "hvy_qpot_includes.h"

int main(int argc, char *argv[])  {
int todo,sm_lev;
int prompt;
Real ssplaq,stplaq;
double dtime;

int i;

initialize_machine(argc,argv);


 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

        dtime = -dclock();
 
	/* Debug */
	/**gaugefix(8,(Real)1.5,500,1e-7,
	   F_OFFSET(staple),F_OFFSET(diag),0,NULL,NULL,0,NULL,NULL);**/

	/* fix to axial gauge */
	if( startflag != CONTINUE){
	    ax_gauge();
	    tot_smear = 0;
	}

	/* Compute unsmeared simple, i.e one-plaquette, glueball operators */
	if( no_smear_level > 0 ){
	    gball_simp(tot_smear);
	}

	/* Loop over the different smearing levels */
	for(sm_lev=0; sm_lev < no_smear_level; sm_lev++ ){

	    /* Do the smearing iterations */
	    for(todo=smear_num[sm_lev]; todo > 0; --todo ){
		smearing();
	    }
	    if(this_node==0)printf("SMEARING COMPLETED\n"); 
	    tot_smear += smear_num[sm_lev];

	    /* Compute simple, i.e one-plaquette, glueball operators */
	    gball_simp(tot_smear);

	    /* Compute on-axis time-like Wilson loops */
	    /** w_loop1(tot_smear); **/

	    /* Compute on-axis time-like hybrid potential loops 
	     and on-axis time-like Wilson loops */
	    hybrid_loop1(tot_smear);

	    /* Compute off-axis time-like Wilson loops, if desired */
	    if( off_axis_flag == 1 ){
		w_loop2(tot_smear);
	    }
	}
        if(this_node==0)printf("RUNNING COMPLETED\n");

        dtime += dclock();
        if(this_node==0){
            printf("Time = %e seconds\n",dtime);
        }
        fflush(stdout);
	dtime = -dclock();

	/* save lattice if requested */
	if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
	}
    }
    return 0;
}
