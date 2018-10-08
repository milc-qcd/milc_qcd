/***************** control.c *****************************************/

/* Measure the disconnected current density correlator as a function of r    */

/* MIMD version 7 */

/* 04/03/15 C. DeTar */

#define CONTROL
#include "rcorr_includes.h"

#ifdef TIME_TO_TIME
#define rcorr rcorr_time
#define print_result print_result_time
#endif

#ifdef T2TFRMPT2PT
#define rcorr rcorr_t2tfrmpt2pt
#define print_result print_result_time
#endif

int
main(int argc, char *argv[])
{
  int prompt;
  complex **qin_sloppy, **qin_diff; // q for quantity in this case current density
  Real *q, **qblock, **q2block;
  double starttime, endtime;
  
  int jflav, k, j, ib;
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
  
    /* Create qin arrays of current densities with each array associated with a random source*/
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

    /* Create q array of the means of the current densities and 
                array of their variances over blocks of random sources */
    qblock = (Real **)malloc(sizeof(Real *)*param.nblock);
    q2block = (Real **)malloc(sizeof(Real *)*param.nblock);
    for(ib = 0; ib < param.nblock; ib++){
#if defined(TIME_TO_TIME) || defined(T2TFRMPT2PT)
      qblock[ib]=(Real *)malloc(sizeof(Real)*nt);
      q2block[ib]=(Real *)malloc(sizeof(Real)*nt);
      for(int t=0;t<nt;t++){
	qblock[ib][t]=0;
	q2block[ib][t]=0;
      }
#else
      qblock[ib] = create_r_field();
      q2block[ib] = create_r_field();
#endif
    }

    /* Read all the data for flavor "jflav" */
    /* Accumulate the result in qin, weighted by the charge  */
    for(jflav = 0; jflav < param.nflav; jflav++){
      accumulate_current_density(param.fname_sloppy[jflav], qin_sloppy, param.charges[jflav], 
				 &param.mass[jflav], param.nrand_sloppy);
      if(strstr(param.fname_diff[jflav],"none") == NULL)
	accumulate_current_density(param.fname_diff[jflav], qin_diff, param.charges[jflav], 
				   &param.mass[jflav], param.nrand_diff);
    }

#ifdef OPT    
    compute_current_density_and_print(qin_sloppy, param.nrand_sloppy, qin_diff, param.nrand_diff, param.fname_curr);
#endif

    /* Calculate the density-density correlator qblock for each random source blocking size */
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

#if !defined(TIME_TO_TIME) && !defined(T2TFRMPT2PT)
    /* Symmetrize over hypercubic group transformations */
    for(ib = 0; ib < param.nblock; ib++)
      symmetrize(qblock[ib], q2block[ib]);
#endif
    /* Extrapolate to infinite block size and 
       write the results to the specified file */
    print_result(qblock, q2block, param.nblock, param.block_size);

    /* Free qblock array */
    for(ib = 0; ib < param.nblock; ib++){
#if !defined(TIME_TO_TIME) && !defined(T2TFRMPT2PT)
      destroy_r_field(qblock[ib]);
      destroy_r_field(q2block[ib]);
#else
      free(qblock[ib]);
      free(q2block[ib]);
#endif	  
    }

    free(qblock);
    free(q2block);
  
  } /* readin(prompt) */

  node0_printf("RUNNING COMPLETED\n");
  endtime=dclock();
    
  node0_printf("Time = %e seconds\n",(double)(endtime-starttime));

  normal_exit(0);
  return 0;
}

/* control.c */
