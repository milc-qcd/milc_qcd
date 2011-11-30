/************************* control.c *******************************/
/* MIMD version 7 */
/* A utility for fixing the gauge field and doing a translation    */

#define CONTROL
#include "gauge_utilities_includes.h"	/* definitions files and prototypes */
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int prompt;
  double dtime, dclock();
  int overrelax = 1.5;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

  dtime = -dclock();

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* gaugefix if requested */
    if( param.fixflag == COULOMB_GAUGE_FIX){
      gaugefix(TUP, overrelax, 20000, param.gauge_fix_tol);
      if(this_node==0)printf("FIXED TO COULOMB GAUGE\n");
      fflush(stdout);
    }
    else if( param.fixflag == LANDAU_GAUGE_FIX){
      gaugefix(8, overrelax,600, param.gauge_fix_tol);
      if(this_node==0)printf("FIXED TO LANDAU GAUGE\n");
      fflush(stdout);
    }

    /* translate the lattice if requested */
    shift_gauge(param.rshift);

    /* Introduce a boundary twist if requested */
    momentum_twist_site(param.bdry_phase,+1);

    /* save lattice if requested */
    if( param.saveflag != FORGET )
      save_lattice( param.saveflag, param.savefile, param.stringLFN );
    
  }	/* end loop over configurations */
  node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
  
  dtime += dclock();
  if(this_node==0){
    printf("Time = %e seconds\n",dtime);
  }
  fflush(stdout);
  
  normal_exit(0);
  return 0;
}
