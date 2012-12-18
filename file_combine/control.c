/***************** control.c ************************************/

/* Extract extended source from staggered or Wilson propagator 		*/
/* MIMD version 7 */

/* Modifications ...
   
 */

#define CONTROL
#include "file_combine_includes.h"
#include <string.h>

#ifndef HAVE_QIO
REQUIRES QIO FOR ALL FILES TREATED HERE
#endif

int main(int argc, char *argv[])
{
  int prompt;
  double starttime, endtime;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

  starttime=dclock();
    
  /* set up */
  prompt = setup();
  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;

    total_iters=0;
    
    if(this_node==0)printf("END OF HEADER\n");
    
    /* Loop over quarks */

    combine_files(param.nfile, param.file_type,
		  param.ncolor, param.nspin, 
		  param.t0,
		  param.startflag, param.startfile,
		  param.coeff, param.saveflag,
		  param.savetype, param.savefile);
    
    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

  } /* readin(prompt) */

  return 0;
}
