/************************ control.c ******************************/
/* MIMD version 6 */
/* Main procedure for SU3 wilson spectrum, including hybrids */

#define CONTROL
#include "wi_hyb_includes.h"

int main(int argc,char *argv[]){
int prompt;
int m_iters,spect_iters;
int disp_iters;
double dtime;

 initialize_machine(argc,argv);
 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){
	dtime = -dclock();

	if(source_start==0){
            /* generate a pseudofermion configuration */
            boundary_flip(MINUS);
            m_iters = f_measure2();
	    boundary_flip(PLUS);
	}
	gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL,
		 F_OFFSET(mp),F_OFFSET(ttt),
		 0,NULL,NULL,0,NULL,NULL);
	boundary_flip(MINUS);
        spect_iters = spectrum_hybrids();
	boundary_flip(PLUS);


	fflush(stdout);
        if(this_node==0)printf("RUNNING COMPLETED\n");
	if(this_node==0)printf("cg/mr iters for spectrum = %d\n", spect_iters);

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);
    }
    return 0;
}
