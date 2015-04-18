/***************** control.c *****************************************/

/* Measure the disconnected current density correlator as a function of r    */

/* MIMD version 7 */

/* 04/03/15 C. DeTar */

#define CONTROL
#include "rcorr_includes.h"

int 
main(int argc, char *argv[])
{
  int prompt;
  complex **qin;
  Real *q;
  double starttime, endtime;
  
  int jflav, k;
  int key[4] = {1,1,1,1};  /* 4D Fourier transform */
  

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();

  starttime=dclock();

  /* set up */
  prompt = setup();

  /* Set up for Fourier transform in all directions */
  if(prompt != 2)setup_restrict_fourier(key, NULL);

  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;   /* For testing */
  
    /* Create qin array */
    qin = (complex **)malloc(sizeof(complex *)*param.nrand);
    if(qin == NULL){
      node0_printf("main: No room for qin\n");
      terminate(1);
    }
    for(k = 0; k < param.nrand; k++){
      qin[k] = create_c_array_field(NMU);
      if(qin[k] == NULL){
	node0_printf("main: No room for qin[%d]\n",k);
	terminate(1);
      }
    }

    for(jflav = 0; jflav < param.nflav; jflav++){

      /* Allocate space and read all the data for flavor "jflav" */
      /* Accumulate the result in qin, weighted by the charge  */

      accumulate_current_density(param.fname[jflav], qin, param.charges[jflav], 
				 &param.mass[jflav], param.nrand);
    }
    
    /* Calculate the density-density correlator q */
    q = rcorr(qin, param.nrand);

    /* Destroy qin array */
    for(k = 0; k < param.nrand; k++)
      destroy_c_array_field(qin[k], NMU);
    free(qin);

    /* Symmetrize over hypercubic group transformations */
    symmetrize(q);

    /* Write the results to the specified file */
    print_result(q, param.nrand);

    destroy_r_field(q);
  
  } /* readin(prompt) */

  node0_printf("RUNNING COMPLETED\n");
  endtime=dclock();
    
  node0_printf("Time = %e seconds\n",(double)(endtime-starttime));

  normal_exit(0);
  return 0;
}

/* control.c */
