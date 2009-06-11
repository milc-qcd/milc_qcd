/***************** control_ext_src.c ************************************/

/* Extract extended source from staggered or Wilson propagator 		*/
/* MIMD version 7 */

/* Modifications ...
   
 */

#define CONTROL
#include "ext_src_includes.h"
#include <string.h>
#ifdef HAVE_QDP
#include "lattice_qdp.h"
#include <qdp.h>
#endif

#ifndef HAVE_QIO
REQUIRES QIO FOR OUTPUT SOURCE FILE
#endif

int main(int argc, char *argv[])
{
  int prompt;
  int i;
  double starttime, endtime;
  int slice[4] = {0,0,0,0};
  int key[4] = {1,1,1,0};  /* We do all time slices, regardless */
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();
  /* loop over input sets */

  /* Initialize the FT */
  setup_restrict_fourier(key, slice);

  while( readin(prompt) == 0){
    
    starttime=dclock();
    
    total_iters=0;
    
    if(this_node==0)printf("END OF HEADER\n");
    
    /* Loop over quarks */

    for(i=0; i<param.num_qk; i++){

      node0_printf("Quark propagator %d\n",i);
      
      if(param.qk_type[i] == CLOVER_TYPE) /* Input Dirac propagator */
	{
	  extract_wprop_to_w_source(param.startflag_w[i], 
				    param.startfile_w[i], 
				    param.num_t0[i],
				    &param.dst_wqs[i][0],
				    &param.snk_wqs[i],
				    param.snk_gam[i]);
	}
      
      else /* Input KS propagator */
	{
	  if(param.dst_type[i] == CLOVER_TYPE) /* Naive extended source */
	    extract_ksprop_to_w_source(param.startflag_ks[i], 
				       param.startfile_ks[i], 
				       param.num_t0[i],
				       &param.dst_wqs[i][0],
				       &param.snk_wqs[i],
				       param.snk_gam[i]);
	  else /* KS extended source */
	    extract_ksprop_to_ks_source(param.startflag_ks[i], 
					param.startfile_ks[i], 
					param.num_t0[i],
					&param.dst_ksqs[i][0],
					&param.snk_ksqs[i]);
	}
    }
      
    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

  } /* readin(prompt) */

  return 0;
}
