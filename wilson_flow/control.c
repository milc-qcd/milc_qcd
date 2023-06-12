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
    run_gradient_flow();

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
