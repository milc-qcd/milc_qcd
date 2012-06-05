/* ************************************************************	*/
/*								*/
/* 			    control_convert.c	   		*/
/*								*/
/* Main procedure for U(1) gauge field conversion       	*/
/* This routine just reads in and writes gauge fields		*/
/*								*/
/* Last updated on 1/05/12 by S. Gottlieb			*/
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
  /*register int i; 
  register site *s; */

  /* Setting up of lattice & neighbors */
  initialize_machine(&argc,&argv);

  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv)==1) terminate(1);
  g_sync();
  prompt=setup();

  /* Main procedure */
  node0_printf("Starting U(1) gauge config conversion run ...\n");
  while(readin(prompt)==0){

	/* Mark time */
	dtime=-dclock();

	g_sync();

	/* save lattices */
	save_u1_lattice(save_u1flag,save_u1file);

	/* check a few things */
	u1plaq(&splq,&tplq);
	plp=u1ploop();
	node0_printf("\nu1-ploop = ( %e, %e )  u1-(s,t)plaq = ( %e, %e )\n",
		plp.real,plp.imag,splq,tplq);

	/* Mark time */
	dtime += dclock();
	node0_printf("U(1) conversion completed!\n");
	node0_printf("Time = %e seconds\n",dtime);
	fflush(stdout);

  } /* while-ends */

  return(0);

} /* end of main() */

/* ************************************************************	*/
