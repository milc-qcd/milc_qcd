/* ************************************************************	*/
/*								*/
/* 			    U1.C		   		*/
/*								*/
/* Main procedure for U(1) Gauge-fixed field generation       	*/
/*								*/
/* Last updated on 07.11.07					*/
/*								*/
/* ************************************************************	*/

#define CONTROL

#include "include_u1g.h"

int main(int argc,char *argv[])
{

  complex plp;
  double dtime;
  Real splq,tplq;
  int prompt;

  /* Setting up of lattice & neighbors */
  initialize_machine(&argc,&argv);

  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv)==1) terminate(1);
  g_sync();

  /* Mark time */
  dtime=-dclock();

  prompt=setup();

  /* Prepare the U(1) field generator.  The momentum table and the FFT
     setup do not change from one configuration to the next, so this is
     done only once. */
  initialize_u1_gauge_field(GAUGE, SCHEME);
  node0_printf("Momentum table prepared!\n"); fflush(stdout);

  /* Main procedure */
  node0_printf("Starting U(1) gauge config generation run ...\n");
  while(readin(prompt)==0){

	/* Generate the U(1) vector potential.  We own the storage for
	   u1_A; see setup() */

	generate_u1_gauge_field(u1_A);

	/* save lattices */
	save_u1_lattice(save_u1flag,save_u1file);

	/* check a few things */
	u1plaq(&splq,&tplq,echarge);
	plp=u1ploop(echarge);
	node0_printf("\nu1-ploop = ( %e, %e ) with charge %f\n", plp.real,plp.imag, echarge);
	node0_printf("u1-(s,t)plaq = ( %e, %e ) with charge %f\n", splq,tplq, echarge);
	double sLink, tLink;
	u1avlink(&sLink, &tLink, echarge);
	node0_printf("u1-(s,t)link = ( %e, %e ) with charge %f\n", sLink, tLink, echarge);

  } /* while-ends */

  finalize_u1_gauge_field();
  destroy_u1_A_field(u1_A);

  /* Mark time */
  dtime += dclock();
  node0_printf("RUNNING COMPLETED\n");
  node0_printf("Time = %e seconds\n",dtime);
  fflush(stdout);

  return(0);

} /* end of main() */

/* ************************************************************	*/

