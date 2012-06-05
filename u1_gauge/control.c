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
  Real splq,tplq,vnorm;
  int prompt,key[4],slc[4];
  register int i,dir;
  register site *s;

  /* Setting up of lattice & neighbors */
  initialize_machine(&argc,&argv);

  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv)==1) terminate(1);
  g_sync();
  prompt=setup();

  /* Initialize fft */
  key[XUP]=1; key[YUP]=1;
  key[ZUP]=1; key[TUP]=1;
  slc[XUP]=0; slc[YUP]=0;
  slc[ZUP]=0; slc[TUP]=0;
  setup_restrict_fourier(key,slc);
  vnorm=(Real)sqrt((double)volume);

  /* Main procedure */
  node0_printf("Starting U(1) gauge config generation run ...\n");
  while(readin(prompt)==0){

	/* Mark time */
	dtime=-dclock();

	/* Generating gauge config in mom space
	   under coulomb and global gaugefixing */

	u1gf = (complex *)malloc(sites_on_node*4*sizeof(complex));
	momgauge(u1gf);

	/* FFT mom space gauge config to coord space */

	restrict_fourier_field(u1gf, 4*sizeof(complex), BACKWARDS);

	FORALLSITES(i,s){
	  FORALLUPDIR(dir){
	    u1_A[4*i+dir] = u1gf[4*i+dir].real/vnorm;
	  }
	}

	free(u1gf);

	/* save lattices */
	save_u1_lattice(save_u1flag,save_u1file);

	/* check a few things */
	u1plaq(&splq,&tplq);
	plp=u1ploop();
	node0_printf("\nu1-ploop = ( %e, %e )  u1-(s,t)plaq = ( %e, %e )\n",
		plp.real,plp.imag,splq,tplq);

	/* Mark time */
	dtime += dclock();
	node0_printf("U(1) running completed!\n");
	node0_printf("Time = %e seconds\n",dtime);
	fflush(stdout);

  } /* while-ends */

  return(0);

} /* end of main() */

/* ************************************************************	*/

