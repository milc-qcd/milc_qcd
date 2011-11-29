/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/control_form.c,v 1.3 2011/11/29 18:04:40 detar Exp $ **/

/*

   Calculation of heavy to light and heavy to heavy form
   form factors, using HQET quarks for the bottom quarks.

*/



#define CONTROL
#include "hqet_light_includes.h"


#ifdef DEBUGDEF
#include DEBUGDEF
#endif


int main(int argc,char **argv) 
{
  int prompt;
  int dummy ; 
  /*---------- start of the calculation ----------  */
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  g_sync();
  /* set up */
  prompt = setup_hqet_form() ;
  /***  more work, I need to add the prompt in the correct MILC style */
  dummy = readin(prompt) ;
  gamma_initialized = 0;

  setup_control_hqet_form();


  calc_hqet_light_form(); 

  normal_exit(0);

  return 0 ;
} /******** end of the main program **********/





