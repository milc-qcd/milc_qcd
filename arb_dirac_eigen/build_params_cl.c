/******** build_params_cl.c *********/
/* MIMD version 6 */
/* clover action */


#include "arb_dirac_eig_includes.h"

void build_params(Real mass_0)
{
int j;
Real cpauli[4];
Real fpauli[4];
Real boost_term;



/* lambda and rho terms */
	for(j=1;j<5;j++) {
        lambda[j]= 0.0;
        rho[j]= 0.0;
	}
        lambda[0]= 4.0 + mass_0;
        lambda[1]=rho[1]= -0.5;


        for(j=0;j<5;j++)
		if(this_node==0) printf("lambda[%d] = %e\n",j,lambda[j]);
        for(j=1;j<5;j++)
		if(this_node==0) printf("rho[%d] = %e\n",j,rho[j]);

/* to renormalize to MILC conventions */
        clover_term = 0.5;



 /* boost clover term */ 
	boost_term=1.0;
        clover_term *= boost_term;
if(this_node==0)printf("boosting clover term by %e to %e\n",boost_term,
clover_term);

if(this_node==0)printf("clover normalization %e\n",clover_term);


} /* build_params */
