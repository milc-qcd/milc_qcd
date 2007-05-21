/**************************** meson_spectrum.c ****************************/
/* MIMD version 6 */
/*
       >>>>>  light-light spectrum functions <<<<<

    This module calculates the degenerate meson 
    spectrum.


 */

#include "prop_form_includes.h"

/**


   Subroutine arguments
      quark          :: site structure pointer (spin_wilson_vector) to the FULL quark propagator
      t_source       :: the time slice of the quark propagator
      what_to_do     :: flag to tell the routine what to do
      ikap           :: which kappa value is being used
      num_kap        :: the total number of kappa values


   This code uses the GLOBAL site structure variable.
   wilson_propagator quark_store ;



**/

void meson_spectrum(field_offset quark,int t_source, int ikap, int num_kap, int what_to_do,char *filename)
{
  int color ; 
  Real norm_fac[10];
  static double_complex *pmes_prop[MAX_KAPPA][10]; /** storage for meson correlators **/
  static char *mes_kind[10] = {"PION","PS505","PS055","PS0505",
				 "RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};

  int num_prop ;
  Real  space_vol = (Real)(nx*ny*nz);
  int t ;
  double t_total ; 
  static int setup_called = 0 ;
  int i ; 
  FILE *fp;

  /***--****************************************--****/

  if( what_to_do == SETUP_CORR )
  {
    /*** storage for the meson correlators *****/
    for(i=0;i<num_kap;i++)
      for(num_prop=0;num_prop<10;num_prop++)
      {
	pmes_prop[i][num_prop] = 
	  (double_complex *)malloc(nt*sizeof(double_complex));
	for(t=0;t<nt;t++)
	{
	  pmes_prop[i][num_prop][t] = dcmplx(0.0,0.0);
	}
      }
      
    IF_VERBOSE_ON(1)
      printf("The MESON correlators have been setup \n") ;

    setup_called = SETUP_CORR ;

    return ; 
  }   /*** end of set up code ***/


  /** the setup code should have been called to get to this part of the code ***/
  assert( setup_called == SETUP_CORR ) ;
  assert( ikap >= 0 && ikap < num_kap ) ;
 /*
  *   calculate the meson spectrum 
  */


  if( what_to_do == CALCULATE_SPECTRUM  )
  {
    w_meson_site(quark, quark, pmes_prop[ikap]);
    return ; 
  }



  if( what_to_do == WRITE_RESULTS  )
  {

    if(saveflag_HH2_GL == SAVE_ASCII)
      {
	fp = fopen(filename,"w");
	if(fp == NULL)
	  {
	    printf("meson_spectrum: Can't open %s for output\n",filename);
	    terminate(1);
	  }

	for(num_prop=0;num_prop<10;num_prop++) 
	  norm_fac[num_prop] = space_vol;
	
	norm_fac[4] *= 3.0;
	norm_fac[5] *= 3.0;
	norm_fac[8] *= 3.0;
	norm_fac[9] *= 3.0;
	
	
	/** sum the correlators over the nodes ****/
	for(i=0;i<num_kap;i++)
	  for(num_prop=0;num_prop<10;num_prop++)
	    for(t=0; t<nt; t++)
	      {
		g_doublesum( &pmes_prop[i][num_prop][t].real );
		pmes_prop[i][num_prop][t].real  /= norm_fac[num_prop];
		g_doublesum( &pmes_prop[i][num_prop][t].imag );
		pmes_prop[i][num_prop][t].imag  /= norm_fac[num_prop];
		
	      }
	
	if(this_node==0)fprintf(fp,"DEGENERATE_MESON_SPECTRUM START\n") ;
	/* print meson propagators */
	for(i=0;i<num_kap;i++)
	  for(num_prop=0;num_prop<10;num_prop++)
	    for(t=0; t<nt; t++)
	      {
		if(this_node == 0)
		  {
		    fprintf(fp,"POINT%s %d %d  %e %e\n",mes_kind[num_prop],t,i,
 		      (double)pmes_prop[i][num_prop][(t+t_source)%nt].real,
		      (double)pmes_prop[i][num_prop][(t+t_source)%nt].imag);
		  }
	      }
	
	
	if(this_node==0)fprintf(fp,"DEGENERATE_MESON_SPECTRUM END\n") ;
	g_sync() ; 
	fclose(fp);

      }
	

  /**
   **   free up the memory
   **/

    for(i=0;i<num_kap;i++)
      for(num_prop=0;num_prop<10;num_prop++)
      {
	free(pmes_prop[i][num_prop] ) ; 
      }


/**  IF_VERBOSE_ON(2)
    printf("Time to calculate meson spectrum = %g sec\n",dclock() - t_total) ;  **/
    
    return ;
  }  /*** end of the final segment of the code ****/



  if( this_node == 0 )
  {
    printf("meson_spectrum:ERROR what_to_do = %d out of range \n",what_to_do) ; 
    terminate(1) ; 
  }


} /**** end of the light meson spectrum function ***/
