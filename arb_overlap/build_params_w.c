/******** build_params_w.c *********/
/* MIMD code version 7 */
/* Good old Wilson action */
#include "arb_ov_includes.h"

void build_params(Real mass_0)
{
int j;
Real cpauli[4];
Real fpauli[4];
Real mass_00;



/* lambda and rho terms */
	for(j=1;j<5;j++) {
        lambda[j]= 0.0;
        rho[j]= 0.0;
	}
	lambda[0]= 4.0 + mass_0;
	lambda[1]=rho[1]= -0.5;

        for(j=0;j<5;j++) node0_printf("lambda[%d] = %e\n",j,lambda[j]);
        for(j=1;j<5;j++) node0_printf("rho[%d] = %e\n",j,rho[j]);

        cpauli[0]= 0.0;
	fpauli[0]=cpauli[0]/16.0;

/* to renormalize to MILC conventions */
	clover_term= -8.0*fpauli[0];
if(this_node==0)printf("clover normalization %e\n",cpauli[0]);


} /* build_params */
