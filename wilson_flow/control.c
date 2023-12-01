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

#ifdef DEBUG_FIELDS
void dump_double_lattice();
#endif

int
main( int argc, char **argv )
{
  /* control variables */
  int prompt;
  double dtime, dtimec, dclock();
  int i;
#ifdef SPHALERON
  double dtimebulk,dtimebdry; 
  Real q_bulk[3];
#endif

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

    /* integrate the flow */
#ifdef REGIONS
    run_gradient_flow( FULLVOL );
    // run_gradient_flow();
#else
    run_gradient_flow();
#endif
    /* Save lattice if requested */
    if( saveflag != FORGET )
      save_lattice( saveflag, savefile, stringLFN );

#ifdef DEBUG_FIELDS
      dump_double_lattice();
#endif

    /* Stop and print timer for this configuration */
    dtimec += dclock();
    node0_printf("Time to complete flow = %e seconds\n", dtimec);
    fflush(stdout);

#ifdef SPHALERON

#ifdef DEBUG_BLOCKING
    test_blocking();
    normal_exit(0);
#else
#ifndef HALF_LATTICE_TEST

    /* Start timer for bulk flow (doesn't include 4D preflow time) */
    dtimebulk = -dclock();
    bulk_flow( q_bulk );
    /* Stop and print timer for bulk flow  */
    dtimebulk += dclock();
    node0_printf("Time to complete bulk flow = %e seconds\n", dtimebulk);
    fflush(stdout);
#else
    report_bulk( stoptime, q_bulk );
#endif
    /* Start timer for bdry flow (doesn't include 4D pre- and bulk-flow time) */
    dtimebdry = -dclock();
    bdry_flow( q_bulk );
    /* Stop and print timer for this configuration */
    dtimebdry += dclock();
    node0_printf("Time to complete bdry flow = %e seconds\n", dtimebdry);
    fflush(stdout);
#endif

#endif
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
