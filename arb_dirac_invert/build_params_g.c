/******** build_params_g.c *********/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */

#include "arb_dirac_inv_includes.h"

/* Gaussian BoC parameterization */
static Real fit_params[9][3]={
{0.0, 0.0, 0.0},
{-2.673198e+00, -7.092352e-01, -1.491836e-02},
{-3.493561e+00, -1.007043e+00, -4.207975e-02},
{-4.211529e+00, -1.123061e+00, -8.036113e-02},
{-4.873246e+00, -1.180296e+00, -1.169135e-01},
{-1.753177e+00, -9.475685e-01, 3.147423e-03},
{-3.496692e+00, -1.030516e+00, -5.020443e-02},
{-4.762544e+00, -1.016991e+00, -1.047275e-01},
{-5.875998e+00, -9.306595e-01, -1.633547e-01},
};


void build_params(Real mass_0)
{
int j,k;
Real chm,shm;
Real rn,as,d;
Real cpauli[4];
Real fpauli[4];
Real an,ax,rx;
Real c123,c132;

Real rt,d1,d2,c40,c01,numer;
Real p[3];

Real boost_term;
Real fudge,ans;
void fpoly(Real x,Real *p, int j);

fudge=1.0/0.952;
boost_term=1.2;


/* lambda and rho terms */
        for(j=1;j<5;j++) {
        fpoly(mass_0,p,3);

        ans=0.0;
        for(k=0;k<3;k++) ans += p[k]*fit_params[j][k];
        lambda[j]= -fudge*exp(ans);

        ans=0.0;
        for(k=0;k<3;k++) ans += p[k]*fit_params[j+4][k];
        rho[j]= -fudge*exp(ans);


        }

/* lambda[0] from pole */
if(mass_0>=0.0){
if(this_node==0) printf(" computing lambda[0] by root finding\n");
        chm=0.5*(exp(mass_0)+exp(-mass_0));
        shm=0.5*(exp(mass_0)-exp(-mass_0));

        rn=2.0*rho[1]+12.0*rho[2]+24.0*rho[3]+16.0*rho[4];
        as=2.0*lambda[1]+12.0*lambda[2]+24.0*lambda[3]+16.0*lambda[4];
        d=8.0*lambda[1]+24.0*lambda[2]+32.0*lambda[3]+16.0*lambda[4];

        lambda[0]= -rn*shm - (chm-1.0)*as -d;
}
else{
if(this_node==0) printf(" computing lambda[0] by extrapolating\n");
/* this attempts to produce an output bare mass_0 equal to the input one */
an=lambda[0]+8.0*lambda[1]+24.0*lambda[2]+32.0*lambda[3]+16.0*lambda[4];
ax=2.0*lambda[1]+12.0*lambda[2]+24.0*lambda[3]+16.0*lambda[4];
as=6.0*lambda[1]+12.0*lambda[2]+8.0*lambda[3];
rx=2.0*rho[1]+12.0*rho[2]+24.0*rho[3]+16.0*rho[4];
rx= -rx;
lambda[0] = rx*mass_0-as-ax;
}



	for(j=0;j<5;j++) if(this_node==0)printf("lambda[%d] = %e\n",j,lambda[j]);
	for(j=1;j<5;j++) if(this_node==0)printf("rho[%d] = %e\n",j,rho[j]);



/* clover term normailzed so mb=m0 */
if(mass_0>=0.0001){
if(this_node==0) printf(" computing  clover term by root finding\n");
/* clover term from solution of mb=m0 */
        for(j=0;j<4;j++) cpauli[j]=0.0;
 
        an=lambda[0]+8.0*lambda[1]+24.0*lambda[2]+32.0*lambda[3]+16.0*lambda[4];
        rt=4.0*rho[2]+16.0*rho[3]+16.0*rho[4];
        rt= -rt; rn= - rn;
 
        numer= -(an+as*(chm-1.0))*(as*shm - rn*chm);
        d1=rn+rt*(chm-1.0);
        d2=an+as*(chm-1.0);
        cpauli[0]=(d1*d1-numer/mass_0)/d2;
        if(this_node==0)printf("C112= %e\n",cpauli[0]);
}
else{
if(this_node==0) printf(" computing  clover term by linear extrapolation\n");
/* c0 = value of c112 at m = 0.01, c40 is value at m=0.4 */        for(j=0;j<4;j++) cpauli[j]=0.0;
        c40= -0.7923;
        c01= -1.024;
	d1=(c40-c01)/0.39;
	d2= c01 - 0.01*d1;
	cpauli[0]= d2+d1*mass_0;
}

	fpauli[0]=cpauli[0]/16.0;
/* to renormalize to MILC conventions */
	clover_term= -8.0*fpauli[0];
if(this_node==0)printf("clover normalization %e\n",cpauli[0]);
 /* boost clover term 
	clover_term *= boost_term;
if(this_node==0)printf("boosting clover term by %e to %e\n",boost_term,
clover_term);
*/


} /* build_params */
