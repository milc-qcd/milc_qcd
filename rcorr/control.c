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
  complex **qin_sloppy, **qin_diff;
  Real **qblock, **q2block;
  double starttime, endtime;
  
  int jflav, k, ib;
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
  
    /* Create qin arrays */
    qin_sloppy = (complex **)malloc(sizeof(complex *)*param.nrand_sloppy);
    if(qin_sloppy == NULL){
      node0_printf("main: No room for qin\n");
      terminate(1);
    }
    for(k = 0; k < param.nrand_sloppy; k++){
      qin_sloppy[k] = create_c_array_field(NMU);
    }

    qin_diff = (complex **)malloc(sizeof(complex *)*param.nrand_diff);
    if(qin_diff == NULL){
      node0_printf("main: No room for qin\n");
      terminate(1);
    }
    for(k = 0; k < param.nrand_diff; k++){
      qin_diff[k] = create_c_array_field(NMU);
    }

    /* Create q array and array of variances of the mean over blocks
       of random sources */
    qblock = (Real **)malloc(sizeof(Real *)*param.nblock); // nblock = # ways to partition the given random vectors into blocks
    q2block = (Real **)malloc(sizeof(Real *)*param.nblock);// variance over blocks of given size.
    for(ib = 0; ib < param.nblock; ib++){
      qblock[ib] = create_r_field();
      q2block[ib] = create_r_field();
    }

    for(jflav = 0; jflav < param.nflav; jflav++){

      /* Read all the data for flavor "jflav" */
      /* Accumulate the result in qin, weighted by the charge  */

      accumulate_current_density(param.fname_sloppy[jflav], qin_sloppy, param.charges[jflav], 
				 &param.mass[jflav], param.nrand_sloppy);
      accumulate_current_density(param.fname_diff[jflav], qin_diff, param.charges[jflav], 
				 &param.mass[jflav], param.nrand_diff);
    }

    /* Calculate the density-density correlator qblock for each random source blocking size (qblock = array of correlator fields) */
    /* Calculate the variance q2block over random source blocks */
    rcorr(qblock, q2block, qin_sloppy, param.nrand_sloppy, qin_diff, param.nrand_diff, 
	  param.nblock, param.block_size);

    /* Destroy qin arrays */
    for(k = 0; k < param.nrand_sloppy; k++)
      destroy_c_array_field(qin_sloppy[k], NMU);
    free(qin_sloppy);

    for(k = 0; k < param.nrand_diff; k++)
      destroy_c_array_field(qin_diff[k], NMU);
    free(qin_diff);

    /* Symmetrize over hypercubic group transformations */
    for(ib = 0; ib < param.nblock; ib++)
      symmetrize(qblock[ib], q2block[ib]);

    /* Extrapolate to infinite block size and 
       write the results to the specified file */
    print_result(qblock, q2block, param.nblock, param.block_size);

    /* Free qblock array */
    for(ib = 0; ib < param.nblock; ib++){
      destroy_r_field(qblock[ib]);
      destroy_r_field(q2block[ib]);
    }
    free(qblock);
  
  } /* readin(prompt) */

  node0_printf("RUNNING COMPLETED\n");
  endtime=dclock();
    
  node0_printf("Time = %e seconds\n",(double)(endtime-starttime));

  normal_exit(0);
  return 0;
}

/* control.c */
