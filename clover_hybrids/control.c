/***************** control.c *********************************/
/* MIMD version 6 */
/* Main procedure for SU3 wilson spectrum, including hybrids */

#define CONTROL
#include "cl_hyb_includes.h"

int main(int argc,char *argv[]){
  int prompt;

  int m_iters,spect_iters;
  double dtime;
  double g_time ; 
  /***  ---------------------------------------- ****/


 initialize_machine(argc,argv);
 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0)
    {
      dtime = -dclock();

	if(source_start==0)
	{
	  /* generate a pseudofermion configuration */
          if (boundary_flag  == ANTI_PERIODIC_IN_TIME ) boundary_flip(MINUS);
/**	  m_iters = f_measure2();  **/
	  if (boundary_flag  == ANTI_PERIODIC_IN_TIME ) boundary_flip(PLUS);
	}

      if( fixflag == COULOMB_GAUGE_FIX)
	{
	  if(this_node == 0) 
	    printf("Fixing to Coulomb gauge\n");

	  g_time = -dclock();

	  gaugefix(TUP,(Real)1.8,500,GAUGE_FIX_TOL,
		   F_OFFSET(mp),F_OFFSET(chi),
		   0,NULL,NULL,0,NULL,NULL); 



	  g_time += dclock();
	  if(this_node==0)printf("Time to gauge fix = %e sec\n",g_time);
	}
      else
      {
	if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
      }

	if (boundary_flag  == ANTI_PERIODIC_IN_TIME ) boundary_flip(MINUS);
      spect_iters = spectrum_hybrids(); 
	if (boundary_flag  == ANTI_PERIODIC_IN_TIME ) boundary_flip(PLUS);


	fflush(stdout);
        if(this_node==0)
	{
	  printf("RUNNING COMPLETED\n");
	  printf("cg/mr iters for spectrum = %d\n", spect_iters);
	  printf("cg/mr iters f_measure2  = %d\n", m_iters);
	}

	dtime += dclock();
	if(this_node==0)
	{
	  printf("Time = %e seconds\n",dtime);
	  printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);
    }



  return 0 ; 
}  /*** end of main ***/
