/********* build_params_0.166.c ********/
/* MIMD version 7 */

/* parameter building for the selected flat action
   fit for r0=1.6 ''best massless''
   and we want the eigenvalue of D(mass=0) so we just linearly
   shift the lambda[0] term

   MIMD code version 6 */

#include "arb_ov_includes.h"

void build_params(float mass_0)
{
  int j;  
  float cpauli[0];
  float fpauli[0];
  float a0_0,a1_0,a2_0,r1_0,r2_0;
/*  float eta,M_add,shm,chm; */
  float M_add;
  static int called=0;
  static float old_mass=0.;

  /* inititalize the flat action values for zero pole mass */
  if (called==0)node0_printf("TEST PROJECTOR planar action \n");
  called=1;
  a0_0=  3.21895;
  a1_0= -0.166666;  r1_0= 0.166666;
  a2_0= -0.078567;  r2_0= 0.055555;
  M_add=mass_0;

  /* lambda and rho terms */
  for(j=3;j<5;j++) {
    lambda[j]= 0.0;
    rho[j]= 0.0;
  }
  lambda[0]= (a0_0+M_add);
  lambda[1]= a1_0;   rho[1]= -r1_0;
  lambda[2]= a2_0;   rho[2]= -r2_0;

  if (mass_0!=old_mass)
  node0_printf("lambda = %e %e %e \n",(double)lambda[0],(double)lambda[1],(double)lambda[2]);
  /*
  for(j=1;j<5;j++) if(this_node==0) printf("rho[%d] = %e\n",j,(double)rho[j]);
  if(this_node==0) printf("clover term frozen to m=0 value\n");
*/
  
  cpauli[0]= -1.278;
  fpauli[0]=cpauli[0]/16.0;
  
  /* to renormalize to MILC conventions */

  clover_term= -8.0*fpauli[0];
  if (mass_0!=old_mass)
  if(this_node==0)printf("clover normalization %e\n",(double)cpauli[0]);
  old_mass=mass_0;

} /* build_params */
