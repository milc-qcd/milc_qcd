/**************************** setup_control_form.c ****************************/
/* MIMD version 6 */
/*
 *  Simple driver routine to collect together all
 *  the functions required to start 
 *
 * The smearing functions in momentum space are used in the 
 * smearing of the quark propagators at the sink.
 * The smearing functions are required at negative momentum.
 *
 *   smear( -p ) = sum_{x} exp( -i p x )  smear(x)
 *  
 *   This corresponds to choosing +1 in "restrict_fourier"
 */


#include "prop_form_includes.h"

#ifdef DEBUGDEF
#include "debug_form.h"
#endif




void setup_control_form()
{
  int i ;
  int k,spin;
  register site *s; 
  double d_time ; 

  d_time = -1.0 * dclock(); 

  /*** set up the FFT routines  *****/
  setup_timeslice_fft() ;

  /**** inintialize the sink smearing function according to wqs ******/
  w_sink_site(F_OFFSET(heavy_smear_func[0]), &wqs_zonked_heavy[0]  ) ; 
  /** warning this code now assumes that the same smearing functions is 
      used for each quark propagator inverted from the origin. This
      should be explicity checked.
  **/

  /***  FFT the zonked smearing function  ***/
  FORALLSITES(i,s)
    {
      s->heavy_smear_func_mom[0] = s->heavy_smear_func[0]; 
    }

  restrict_fourier_site(F_OFFSET(heavy_smear_func_mom[0]), 
			sizeof(complex), 1);


  for(i=0 ; i < no_p_values ; ++i)
  {
    load_smearing(F_OFFSET(seq_smear_func[ i ] ),  seq_smear_file[ i ]  ) ;

  }  /*** end of the loop over the insertion momentum values ***/


  assert( MAXPMOM < 48 ) ;
  /** The FFT work space should be less than a spin_wilson_vector  **/

  /*** Fourier transform the smearing functions ***/

  restrict_fourier_site(F_OFFSET(seq_smear_func[0]),
			no_p_values*sizeof(complex), 1);

  /*** Zero the primary quark propagators ***/

  FORALLSITES(i,s)
    {
      for(k = 0; k < no_zonked_light; k++)
	for(spin=0; spin<4; spin++)
	  clear_wvec(&s->quark_zonked_light[k].d[spin]);
      for(k = 0; k < no_zonked_heavy; k++)
	for(spin=0; spin<4; spin++)
	  clear_wvec(&s->quark_zonked_heavy[k].d[spin]);
      for(spin=0; spin<4; spin++)
	clear_wvec(&s->quark_spectator.d[spin]);
    }


  d_time += dclock(); 
  IF_VERBOSE_ON(1)
    printf("Time in setup_control_form = %g sec\n",d_time) ;

}
