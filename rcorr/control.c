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
  complex *qin[MAXRAND];
  Real *q;
  
  int jrand, nrand = 0;
  int jflav;
  int key[4] = {1,1,1,1};  /* 4D Fourier transform */
  

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  /* We keep a separate field for each random source */
  for(jrand = 0; jrand < MAXRAND; jrand++)
    qin[jrand] = NULL;
  
  /* set up */
  prompt = setup();

  /* Set up for Fourier transform in all directions */
  if(prompt != 2)setup_restrict_fourier(key, NULL);

  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;   /* For testing */
  
    for(jflav = 0; jflav < param.nflav; jflav++){

      /* Allocate space and read all the data for flavor "jflav" */
      /* Accumulate the result in qin, weighted by the charge  */

      accumulate_current_density(param.fname[jflav], qin, param.charges[jflav], 
				 &param.mass[jflav], &nrand);
    }
    
    if(nrand < 2){
      fprintf(stderr, "You need more than 1 random source to compute correlations\n");
      break;
    }

    param.nrand = nrand;
    
    /* Calculate the density-density correlator q */
    q = rcorr(qin, nrand);

    /* Symmetrize over hypercubic group transformations */
    symmetrize(q);

    /* Write the results to the specified file */
    print_result(q, nrand);

    destroy_r_field(q);
  
  } /* readin(prompt) */

  normal_exit(0);
  return 0;
}

/* control.c */
