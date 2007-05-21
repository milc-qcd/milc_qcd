/**************************** setup_control_hqet_form.c ****************************/
/* MIMD version 6 */
/*
 *  Simple driver routine to collect together all
 *  the functions required to start the hqet heavy to light 
 *  form factor code.
 *
 */


#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


void setup_control_hqet_form(void)
{
  int ivel, i;
  int sgn  ;
  double t_total ; 
  register site *s; 
  /*******------------------------------**********/

  t_total = dclock() ; 

  setup_timeslice_fft() ;

  /*** load in all the required smearing functions ****/

  for(ivel=0 ; ivel < novel ; ++ivel)
  {
    load_smearing(F_OFFSET(seq_smear_func[ ivel ] ),  hqet_smear_file[ ivel ]  ) ;

    /*** copy the smearing functions ***/
    FORALLSITES(i,s) 
    {
      s->seq_smear_func_fft[ ivel ] = s->seq_smear_func[ ivel ] ;
    }

  }  /*** end of the loop over the velocites ***/


  /*** Fourier transform the smearing functions ***/
  /**

    f(-p) = \sum_{\x} exp( ipx )  f(x) corresponds to isign = -1 

  **/

  sgn = -1  ;
  restrict_fourier_site(F_OFFSET(seq_smear_func_fft[ 0 ]),
			novel*sizeof(complex), sgn);



  IF_VERBOSE_ON(1)
    printf("setup_control_hqet_form::Total time in this function = %g sec\n",dclock() - t_total );


}   /*** end of setup_control_hqet_form() *****/





