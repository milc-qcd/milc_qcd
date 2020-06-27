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
#if (MILC_PRECISION==1)
    node0_printf("WFLOW: %g %g %g %g\n", 0.0, Et, Es, charge);
#else
    node0_printf("WFLOW: %g %.16g %.16g %.16g\n", 0.0, Et, Es, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
      node0_printf("ADAPT: %g %g %g %g\n", 0.0, stepsize, 0.0, 0.0 );
#else
      node0_printf("ADAPT: %g %.16g %.16g %.16g\n", 0.0, stepsize, 0.0, 0.0 );
#endif
#endif
    fflush(stdout);

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    steps_rejected = 0; // count rejected steps in adaptive schemes
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    is_first_step = 1; // need to know the first step for FSAL
    // set the permutation array for FSAL
    indK[0] = 0; indK[1] = 1; indK[2] = 2; indK[3] = 3;
#endif
    is_final_step = 0;
    flowtime = 0;
    i = 0;
    /* Loop over the flow time */
    while( stoptime==AUTO_STOPTIME || ( flowtime<stoptime && is_final_step==0 ) ) {
//    for( flowtime=stepsize, i=0;
//         stoptime==AUTO_STOPTIME || flowtime<stoptime;
//         flowtime+=stepsize, i++ ) {
      /* Adjust last time step to fit exactly stoptime */
      if( stepsize>stoptime-flowtime && stoptime!=AUTO_STOPTIME ) {
        stepsize = stoptime-flowtime;
        is_final_step = 1;
      }
//      printf("%g\n", stepsize);

      /* Perform one flow step (most of the computation is here) */
      flow_step();

      flowtime += stepsize;
      i++;

      /* Calculate and print current flow output */
      fmunu_fmunu(&Et, &Es, &charge);
#if (MILC_PRECISION==1)
      node0_printf("WFLOW: %g %g %g %g\n", flowtime, Et, Es, charge);
#else
      node0_printf("WFLOW: %g %.16g %.16g %.16g\n", flowtime, Et, Es, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
      node0_printf("ADAPT: %g %g %g %g\n", flowtime, stepsize, dist, local_tol/dist );
#else
      node0_printf("ADAPT: %g %.16g %.16g %.16g\n", flowtime, stepsize, dist, local_tol/dist );
#endif
#endif
      fflush(stdout);

      /* Automatic determination of stoptime:                         */
      /*  t^2 E > 0.45 and d/dt { t^2 E } > 0.35                      */
      /*  Bounds need to be adjusted with scale determination cutoff  */
      if( stoptime==AUTO_STOPTIME ) {

        old_value = new_value;
        new_value = flowtime*flowtime*(Et+Es);
        der_value = flowtime*(new_value-old_value)/stepsize;

        if( new_value > 0.45 && der_value > 0.35 )
          break;
      } /* end: auto stoptime */

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
      if( is_final_step==0 ) {
        // adjust step size for the next step except if it is final
        stepsize = stepsize * SAFETY * pow( local_tol/dist, 1/3. );
      }
#endif

    } /* end: flowtime loop */

    /* Save and print the number of steps */
    total_steps = i;
    node0_printf("Number of steps = %i\n", total_steps);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    node0_printf("Number of rejected steps = %i\n", steps_rejected);
#endif
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
