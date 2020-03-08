/************************* control.c *******************************/
/* MIMD version 7 */
/* A utility for fixing the gauge field and doing a translation    */

#define CONTROL
#include "gauge_replicate_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int prompt;
  double dtime, starttime, endtime;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

  starttime = dclock();

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    if(prompt == 2)continue;

    /* Replicate the gauge field */
    if( param.saveflag != FORGET )
      save_replicated_lattice(param.reps, 
			      param.saveflag, param.savefile, 
			      param.stringLFN,
			      QIO_SINGLEFILE, QIO_SERIAL,
			      QIO_ILDGNO, NULL);
    
  }	/* end loop over configurations */

  node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
  
  endtime = dclock();
  if(this_node==0){
    printf("Time = %e seconds\n",(double)(endtime-starttime));
  }
  fflush(stdout);
  
  normal_exit(0);
  return 0;
}
