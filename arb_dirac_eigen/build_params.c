/********* build_params.c ********/
/* MIMD version 6 */

/* parameter building for the selected flat action
   fit for r0=1.6 ''best massless''
   and we want the eigenvalue of D(mass=0) so we just linearly
   shift the lambda[0] term

   MIMD version 6 */

#include "arb_dirac_eig_includes.h"

void build_params(Real mass_0)
{
  int j,k;  
  Real cpauli[4];
  Real fpauli[4];
  Real rn,as,d;
  Real an,ax,rx;
  Real c123,c132;
  Real rt,d1,d2,c40,c01,numer;
  Real a0_0,a1_0,a2_0,r1_0,r2_0;
  Real Rx_0,Ax_0;
  Real eta,M_add,shm,chm;

  /* inititalize the flat action values for zero pole mass */
  a0_0=  2.83444;
  a1_0= -0.17016;  r1_0= 0.17693;
  a2_0= -0.06138;  r2_0= 0.05384;
  Ax_0=2*a1_0+12*a2_0;
  Rx_0=2*r1_0+12*r2_0;
  /* compute residue and mass additive term */
/*
  shm=sinh(mass_0); chm=cosh(mass_0);
  eta=Rx_0*chm-Ax_0*shm;
  M_add=Rx_0*shm-Ax_0*(chm-1);
*/
M_add=mass_0;

  /* lambda and rho terms */
  for(j=3;j<5;j++) {
    lambda[j]= 0.0;
    rho[j]= 0.0;
  }
  lambda[0]= (a0_0+M_add);
  lambda[1]= a1_0;   rho[1]= -r1_0;
  lambda[2]= a2_0;   rho[2]= -r2_0;

  for(j=0;j<5;j++) printf("lambda[%d] = %e\n",j,lambda[j]);
  for(j=1;j<5;j++) printf("rho[%d] = %e\n",j,rho[j]);

  /* clover term normailzed so mb=m0 */
  if(mass_0>=0.0001){
    if(this_node==0) printf(" computing  clover term by root finding\n");
    /* clover term from solution of mb=m0 */
    for(j=0;j<4;j++) cpauli[j]=0.0;
 
    an=lambda[0]+8.0*lambda[1]+24.0*lambda[2]+32.0*lambda[3]+16.0*lambda[4];
    as=2.0*lambda[1]+12.0*lambda[2]+24.0*lambda[3]+16.0*lambda[4];
    rn=2.0*rho[1]+12.0*rho[2]+24.0*rho[3]+16.0*rho[4];
    rt=4.0*rho[2]+16.0*rho[3]+16.0*rho[4];
   rt= -rt; rn= - rn;   /* WATCH FOR POSSIBLE SIGN CONVENTION PROBLEMS */
 
    numer= -(an+as*(chm-1.0))*(as*shm - rn*chm);
    d1=rn+rt*(chm-1.0);
    d2=an+as*(chm-1.0);
    cpauli[0]=(d1*d1-numer/mass_0)/d2;
    if(this_node==0)printf("C112= %e\n",cpauli[0]);
  }
  else{
    if(this_node==0) printf(" computing  clover term by linear extrapolation\n");
    /* c0 = value of c112 at m = 0.01, c40 is value at m=0.4 */        
    for(j=0;j<4;j++) cpauli[j]=0.0;
    c40= -0.7923;
    c01= -1.024;
    d1=(c40-c01)/0.39;
    d2= c01 - 0.01*d1;
    cpauli[0]= d2+d1*mass_0;
  }
printf("clover term frozen to m=0 value\n");
cpauli[0]= -1.029941;

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
