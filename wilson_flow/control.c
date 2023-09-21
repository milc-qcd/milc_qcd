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
  double dtimeb; 
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
#ifdef SPHALERON
    run_gradient_flow_region( FULLVOL );
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
    /* Start timer for this configuration (doesn't include load time) */
    dtimeb = -dclock();
#ifdef DEBUG_BLOCKING
    test_blocking();
    normal_exit(0);
#else
    prepare_bulk_links();
    report_bulk( stoptime );
    /* integrate the flow */
    run_gradient_flow_region( BULK );
    report_bulk( stoptime + stoptime_bulk );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + stoptime_bdry );
    spatial_blocking();
    report_bulk( stoptime + stoptime_bulk + stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 2. * stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 3. * stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 4. * stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 5. * stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 6. * stoptime_bdry );
    run_gradient_flow_region( BOUNDARY );
    report_bulk( stoptime + stoptime_bulk + 7. * stoptime_bdry );
#endif
    /* Stop and print timer for this configuration */
    dtimeb += dclock();
    node0_printf("Time to complete bulk flow = %e seconds\n", dtimeb);
    fflush(stdout);
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
