/**************************** spectrum_light_light.c ****************************/
/* MIMD version 7 */
/*
       >>>>>  light-light spectrum functions <<<<<

    This module calculates the degenerate meson and baryon
    spectrum.


 */

#include "cl_hyb_includes.h"

/**


   Subroutine arguments
      quark          :: site structure pointer (wilson_matrix) to the FULL quark propagator
      t_source       :: the time slice of the quark propagator


   This code uses the GLOBAL site structure variable.
   wilson_propagator quark_store ;



**/

void light_meson_spectrum(int t_source)
{
  int color ; 
  Real norm_fac[10];
  complex *pmes_prop[10]; /** storage for meson correlators **/
  static char *mes_kind[10] = {"PION","PS505","PS055","PS0505",
				 "RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};

  int num_prop ;
  Real  space_vol = (Real)(nx*ny*nz);
  int t ;
  double t_total ; 
  /***--****************************************--****/
  t_total = dclock() ; 



  /*** storage for the meson correlators *****/
  for(num_prop=0;num_prop<10;num_prop++)
  {
    pmes_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
    for(t=0;t<nt;t++)
    {
      pmes_prop[num_prop][t] = cmplx(0.0,0.0);
    }
  }


 /*
  *   calculate the meson spectrum 
  */


    for(color = 0 ; color < 3 ; ++color)
    {

      w_meson_site(F_OFFSET(quark_store.c[color]),
	      F_OFFSET(quark_store.c[color]),pmes_prop);

    } /*** end the loop over colour ***/


    for(num_prop=0;num_prop<10;num_prop++) 
	norm_fac[num_prop] = space_vol;

    norm_fac[4] *= 3.0;
    norm_fac[5] *= 3.0;
    norm_fac[8] *= 3.0;
    norm_fac[9] *= 3.0;


    /** sum the correlators over the nodes ****/
    for(num_prop=0;num_prop<10;num_prop++)
      for(t=0; t<nt; t++)
      {
	g_floatsum( &pmes_prop[num_prop][t].real );
	pmes_prop[num_prop][t].real  /= norm_fac[num_prop];
	g_floatsum( &pmes_prop[num_prop][t].imag );
	pmes_prop[num_prop][t].imag  /= norm_fac[num_prop];
	
      }

    /* print meson propagators */
    for(num_prop=0;num_prop<10;num_prop++)
      for(t=0; t<nt; t++)
      {
	if(this_node == 0)
	{
	  printf("POINT%s %d  %e %e\n",mes_kind[num_prop],t,
		 (double)pmes_prop[num_prop][(t+t_source)%nt].real,
	  (double)pmes_prop[num_prop][(t+t_source)%nt].imag);
	}
      }







  /**
   **   free up the memory
   **/

  for(num_prop=0;num_prop<10;num_prop++)
  {
    free(pmes_prop[num_prop] ) ; 
  }


  IF_VERBOSE_ON(2)
    printf("Time to calculate meson spectrum = %g sec\n",dclock() - t_total) ; 

} /**** end of the light meson spectrum function ***/


/**  Calculate the baryon spectrum and print it to the 
 **  screen
 **
 **  Subroutine arguments
 **   quark     :: site structure pointer (wilson_propagator) to the FULL quark propagator
 **   t_source  :: the time slice of the quark propagator
 **
 **/


void light_baryon_spectrum(field_offset quark, int t_source)
{
  static char *bar_kind[4] = {"PROTON","PROTON0","DELTA","DELTA0"};
  wilson_prop_field *wp;

  int num_prop ;
  Real  space_vol = (Real)(nx*ny*nz);
  int t ;
  int c,i;
  site *s;

  complex *bar_prop[4];
  double t_total ; 
  /***--****************************************--****/

  t_total = dclock() ; 

  for(num_prop=0;num_prop<4;num_prop++)
  {
    if( (bar_prop[num_prop] = (complex *)malloc(nt*sizeof(complex)) ) == NULL )
    {
      if( this_node == 0 ) 
	printf("light_baryon_spectrum:: Error could not reserve the memory for the baryon correlators\n") ; 
      terminate(1); 
    }

    for(t=0;t<nt;t++)
      bar_prop[num_prop][t] = cmplx(0.0,0.0); 
  }

  /*** calculate the baryon spectrum ****/

  /* Create the appropriate field for w_baryon and copy */
  wp = create_wp_field(3);
  for(c = 0; c < 3; c++){
    FORALLSITES(i,s){
      wp->swv[c][i] = ((wilson_propagator *)F_PT(s,quark))->c[c];
    }
  }

  w_baryon(wp, wp, wp, bar_prop );

  destroy_wp_field(wp);

  /**
      print baryon propagators 
   **/

  for(num_prop=0;num_prop<4;num_prop++)
    for(t=0; t<nt; t++)
    {
      g_floatsum( &bar_prop[num_prop][t].real );
      bar_prop[num_prop][t].real  /= space_vol;
      g_floatsum( &bar_prop[num_prop][t].imag );
      bar_prop[num_prop][t].imag  /= space_vol;
    }

  g_sync() ; 


  if(this_node == 0)
  {
    for(num_prop=0;num_prop<4;num_prop++)
      for(t=0; t<nt; t++)
      {
	
	printf("POINT%s %d   %e %e\n",bar_kind[num_prop],t,
	       (double)bar_prop[num_prop][(t+t_source)%nt].real,
	       (double)bar_prop[num_prop][(t+t_source)%nt].imag);
      }
  } /*** end of node 0 ****/


  /**
   **   free up the memory
   **/

      for(num_prop=0;num_prop<4;num_prop++)
      {
	free(bar_prop[num_prop]);
      }

  IF_VERBOSE_ON(2)
    printf("Time to calculate baryon spectrum = %g sec\n",dclock() - t_total) ; 


} /**** end of the light baryon spectrum function ***/
