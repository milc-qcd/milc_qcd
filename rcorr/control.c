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
  complex *qin_sloppy, *qin_diff;
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
  
//    /* Create qin array */
//    qin = (complex **)malloc(sizeof(complex *)*param.nrand);
//    if(qin == NULL){
//      node0_printf("main: No room for qin\n");
//      terminate(1);
//    }
//    for(k = 0; k < param.nrand; k++){
//      qin[k] = create_c_array_field(NMU);
//      if(qin[k] == NULL){
//	node0_printf("main: No room for qin[%d]\n",k);
//	terminate(1);
//      }
//    }

    qin_sloppy = create_c_array_field(NMU);
    qin_diff = create_c_array_field(NMU);

    for(jflav = 0; jflav < param.nflav; jflav++){

      /* Read all the data for flavor "jflav" */
      /* Accumulate the result in qin, weighted by the charge  */

      accumulate_current_density(param.fname_sloppy[jflav], qin_sloppy, param.charges[jflav], 
				 &param.mass[jflav], param.nrand_sloppy);
      accumulate_current_density(param.fname_diff[jflav], qin_diff, param.charges[jflav], 
				 &param.mass[jflav], param.nrand_diff);
    }

    
    /* Calculate the density-density correlator q */
    q = rcorr(qin_sloppy, param.nrand_sloppy, qin_diff, param.nrand_diff);

    /* Destroy qin array */
//    for(k = 0; k < param.nrand; k++)
//      destroy_c_array_field(qin[k], NMU);
//    free(qin);

    destroy_c_array_field(qin_diff, NMU);
    destroy_c_array_field(qin_sloppy, NMU);

    /* Symmetrize over hypercubic group transformations */
    symmetrize(q);

    /* Write the results to the specified file */
    print_result(q);

    destroy_r_field(q);
  
  } /* readin(prompt) */

  node0_printf("RUNNING COMPLETED\n");
  endtime=dclock();
    
  node0_printf("Time = %e seconds\n",(double)(endtime-starttime));

  normal_exit(0);
  return 0;
}

/* control.c */
