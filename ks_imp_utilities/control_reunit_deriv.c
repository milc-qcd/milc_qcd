/************************* control.c *******************************/
/* MIMD version 7 */

/* This code carries out one step in the fermion force calculation,
   namely, the contribution from the derivative of the reunitarized
   link with respect to the preunitarized link and either generates
   the fiducial result or compares with a fiducial result. */

#define CONTROL
#include "ks_imp_utilities_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#endif
#include "params.h"
#ifdef HAVE_GRID
#include "../include/generic_grid.h"
#endif

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
  int prompt;

#ifdef PRTIME
  double dtime;
#endif
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O if needed */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

  double starttime=dclock();
    
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  /* loop over input sets */
  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;
    
    node0_printf("BEGIN\n");
    
    check_reunitarization_derivative( param.ansfile[0], param.ansflag[0],
				      param.ansfile[1], param.ansflag[1] );
    
    node0_printf("RUNNING COMPLETED\n");
    double endtime=dclock();
  
    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    starttime = endtime; /* In case we continue looping over readin */
  
  } /* readin(prompt) */

#ifdef HAVE_GRID
  finalize_grid();
#endif

  normal_exit(0);
  return 0;
}

