/************************* control.c *******************************/
/* Integrates the gauge fields with Wilson or Symanzik flow        */

/* Flag for control file of application */
#define CONTROL

/* Definitions, files, and prototypes */
#include "wilson_flow_includes.h"

/* Additional includes */
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

int
main( int argc, char **argv )
{
  /* control variables */
  int prompt;
  double dtime, dtimec, dclock();

  /* RK integration variables */
  int i;
  double flowtime;

  /* Wilson flow output variables */
  double Et, Es, charge;
  double old_value=0, new_value=0;
  double der_value=0;

  /* Initialization */
  initialize_machine(&argc, &argv);
  if( remap_stdio_from_args(argc, argv) == 1 )
    terminate(1);
  g_sync();

  /* Start application timer */
  dtime = -dclock();

  /* Setup lattice parameters */
  prompt = setup();

  /* Allocate memory for temporary gather matricies */
  for( i=0; i<N_TEMPORARY; i++ )
    tempmat[i] = (su3_matrix *)malloc(sites_on_node * sizeof(su3_matrix));

  /* Loop over configurations */
  while( readin(prompt) == 0 ) {

    /* Start timer for this configuration (doesn't include load time) */
    dtimec = -dclock();

    /* Print flow output column labels */
    node0_printf("#LABEL time Et Es charge\n");
    fflush(stdout);

    /* Calculate and print initial flow output */
    fmunu_fmunu(&Et, &Es, &charge);
    node0_printf("WFLOW %g %g %g %g\n", 0.0, Et, Es, charge);
    fflush(stdout);

    /* Loop over the flow time */
    for( flowtime=stepsize, i=0; 
         stoptime==AUTO_STOPTIME || flowtime<=stoptime; 
         flowtime+=stepsize, i++ ) {

      /* Perform one flow step (most of the computation is here) */
      stout_step_rk();

      /* Calculate and print current flow output */
      fmunu_fmunu(&Et, &Es, &charge);
      node0_printf("WFLOW %g %g %g %g\n", flowtime, Et, Es, charge);
      fflush(stdout);

      /* Automatic determination of stoptime:                         */
      /*  t^2 E > 0.45 and d/dt { t^2 E } > 0.35                      */
      /*  Bounds need to be adjusted with scale determination cutoff */
      if( stoptime==AUTO_STOPTIME ) {

        old_value = new_value;
        new_value = flowtime*flowtime*(Et+Es); 
        der_value = flowtime*(new_value-old_value)/stepsize;

        if( new_value > 0.45 && der_value > 0.35 ) 
          break;
      } /* end: auto stoptime */
    } /* end: flowtime loop */

    /* Save and print the number of steps */
    total_steps = i;
    node0_printf("Number of steps = %i\n", total_steps);
    fflush(stdout);

    /* Save lattice if requested */
    if( saveflag != FORGET )
      save_lattice( saveflag, savefile, stringLFN );

    /* Stop and print timer for this configuration */
    dtimec += dclock();
    node0_printf("Time to complete flow = %e seconds\n", dtimec);
    fflush(stdout);
    
  }/* end: loop over configurations */

  /* Notify user application is done */
  node0_printf("RUNNING COMPLETED\n"); 
  fflush(stdout);

  /* Stop and print application timer */
  dtime += dclock();
  node0_printf("Time = %e seconds\n", dtime); 
  fflush(stdout);

  normal_exit(0);
  return 0;
}
