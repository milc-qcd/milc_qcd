/************************* control.c *******************************/
/* MIMD version 6 */
/* Main procedure for SU3 with dynamical fermions 			*/
/* This application is for measurements on stored lattices.  All the
   machinery for generating lattices is stripped out. */

#define CONTROL
#include "ks_hyb_includes.h"

int main(int argc,char **argv) {
int meascount,todo;
int prompt;
double ssplaq,stplaq;
Real rpbpp,rpbpm;
Real f_energy,f_pressure;	/* fermionic energy and pressure */
Real f_action;			/* fermionic action */
Real rsq;
int m_iters,s_iters,spect_iters,avm_iters,avs_iters;
double dtime;
register int i,dir;
register site *s;

 initialize_machine(argc,argv);
 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

       /* Check input plaquette and unitarity */
        d_plaquette(&ssplaq,&stplaq);
        if(this_node==0){
            printf("CHECK PLAQ: %e %e\n",-ssplaq,-stplaq);fflush(stdout);
        }
        check_unitarity();
        if(this_node==0)printf("Unitarity checked\n"); fflush(stdout);
	dtime = -dclock();

        /* gaugefix and smear links, copy result into original links */
        dtime = -dclock();
        rephase( OFF );
        gaugefix(TUP,(Real)1.8,300,(Real)GAUGE_FIX_TOL,
	   F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),0,NULL,NULL,0,NULL,NULL);
        smear_links( F_OFFSET(link[0]), F_OFFSET(smearlink[0]), staple_weight);
	FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
	    s->link[dir] = s->smearlink[dir];
	}
        rephase( ON );
        /* Check plaquette on fat link lattice */
        d_plaquette(&ssplaq,&stplaq);
        if(this_node==0){
            printf("FAT PLAQ: %e %e\n",-ssplaq,-stplaq);fflush(stdout);
        }

	/* measure the spectrum */
        spect_iters = spectrum_hybrids();
	    
	if(this_node==0)printf("cg iters for spectrum = %d\n", spect_iters);
	fflush(stdout);

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);

	/* save lattice if requested */
        if( saveflag != FORGET ){
          rephase( OFF );
          save_lattice( saveflag, savefile );
          rephase( ON );
        }
    }
    return 0;
}
