/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/propagating_form_factor/control_form.c,v 1.3 2011/11/29 18:49:04 detar Exp $ **/

/*

   Calculation of heavy to light and heavy to heavy form
   form factors, using propagating Wilson quarks.

*/



#define CONTROL
#include "prop_form_includes.h"


#ifdef DEBUGDEF
#include "debug_form.h"
#endif

int main(int argc,char **argv) 
{
  int prompt;
  int dummy ; 
  double total_time ; 
  /*---------- start of the calculation ----------  */
  total_time = dclock() ; 

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  /* set up */
  prompt = setup_h();

  gamma_initialized = 0;
  while( readin(prompt) == 0 )
  {
    setup_control_form();
    calc_heavy_light_form() ;
    IF_MASTER 
    {
      printf("RUNNING COMPLETED\n"); 
      printf("Total time = %g sec\n",  dclock()  -   total_time ); 
    }
  }


  return 0 ;
} /******** end of the main program **********/





