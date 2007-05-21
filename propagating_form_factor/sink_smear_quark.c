/**************************** sink_smear_quark.c ****************************/
/* MIMD version 6 */
/*
   Various routines to smear the light quark
   propagator at the sink.

   Different routines take the quark propagator in 
   different formats.

 */



#include "prop_form_includes.h"
#include <assert.h>

#ifdef DEBUGDEF
#include "debug_form.h"
#endif



/*
 *   Smear the light quark propagator at the sink
 *
 *   quark_smear(x) = \sum_{r) quark(x+r) * smear(r)
 *
 *                  = 1/V \sum_{p) exp( -i p x)   quark(p) * smear(-p)
 *
 *   This function assumes that the quark propagator and smearing functions
 *   are in momentum space. This routine just ties the propagator and smearing
 *   functions together and FFT's back into position space.
 *
 *    The code now includes the volume factors
 *
 *
 *   Subroutine arguments 
 *
 *       quark        :: wilson_vector pointer to the site structure, 
 *                     the quark propagator for a specific
 *                     source colour and spin in momentum space.
 *                     The smeared quark propagator is returned in this 
 *                     data structure.
 *
 *       fft_work_a  :: site structure wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *       fft_work_b  :: site structure wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *      smear_func_mom :: site structure wilson_vector pointer 
 *                         smearing function in momentum space                
 *
 *
 */

void sink_smear_light_quark(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func_mom )
{
  register site *s;
  int i ; 
  double tt ;
  complex nrm ; 

  nrm.real = 1.0 / (nx*ny*nz)   ; nrm.imag = 0.0 ; 

  tt = -1.0*dclock() ; 

  /** multiply the quark propagator by the smearing function in momentum space ****/
  FORALLSITES(i,s)
  {
    c_scalar_self_mult_wvec((wilson_vector *) F_PT(s,quark ),
			     *((complex *) F_PT(s, smear_func_mom)) );

  }


  /** FFT back to real space **/
  restrict_fourier_site(quark, sizeof(wilson_vector), 1); 

  /*** The 1 in the restrict_fourier routine corresponds to the 
       phase factor of exp( -i p x )
   ***/

  tt += dclock() ; 


  IF_VERBOSE_ON(1)
    printf("Time in sink_smear_light_quark = %e sec\n",tt) ;

}




/*
 *  Multiply a wilson vector with a complex number
 *
 *  ans  ----> ans * s 
 *
 *
 */


void c_scalar_self_mult_wvec(wilson_vector *ans, complex s) 
{
  int i , j ; 
  complex z ; 

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 3 ; ++j)
    {
      z = cmul(&(ans->d[ i ].c[ j ]  ),&s ) ;
      ans->d[ i ].c[ j ]  = z ; 
    }


}




void c_scalar_self_mult_spin_wvec(spin_wilson_vector *ans, complex s) 
{
  int i,j,k ; 
  complex z ; 

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 3 ; ++j)
      for(k = 0 ; k < 4 ; ++k)
      {
	z = cmul(&(ans->d[k].d[ i ].c[ j ]  ),&s ) ;
	ans->d[k].d[ i ].c[ j ]  = z ; 
      }


}


/*    >>>>>  wilson_vector data structures <<<<<
 *
 *   Smear the light quark propagator at the sink
 *
 *   quark_smear(x) = \sum_{r) quark(x+r) * smear(r)
 *
 *                  = 1/V \sum_{p) exp( -i p x)   quark(p) * smear(-p)
 *
 *   This function assumes that the quark propagator is in position space
 *   and that the smearing function is already in momentum space.
 *
 *    The code now includes the volume factors
 *
 *
 *   Subroutine arguments 
 *
 *       quark        :: wilson_vector pointer to the site structure, 
 *                     the quark propagator for a specific
 *                     source colour and spin in POSITION space.
 *                     The smeared quark propagator is returned in this 
 *                     data structure.
 *
 *       fft_work_a  :: site structure wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *       fft_work_b  :: site structure wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *      smear_func_mom :: site structure wilson_vector pointer 
 *                         smearing function in momentum space                
 *
 *   Any volume factors are multiplied into the correlators later in the code.
 *
 */


void sink_smear_light_pos_quark(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func_mom )
{
  register site *s;
  int i ; 
  double tt ;
  complex nrm ; 

  nrm.real = 1.0 / (nx*ny*nz)   ; nrm.imag = 0.0 ; 



  tt = -1.0*dclock() ; 

  /*** fourier transform the quark propagator 
    quark( p ) = sum_x exp( i p .x ) quark(x)

    This corresponds to the negative sign in restrict_fourier.
   *****/

  restrict_fourier_site(quark, sizeof(wilson_vector), -1); 

  /** multiply the quark propagator by the smearing function in momentum space ****/
  FORALLSITES(i,s)
  {
    c_scalar_self_mult_wvec((wilson_vector *) F_PT(s,quark ),
			     *((complex *) F_PT(s, smear_func_mom)) );


    c_scalar_self_mult_wvec((wilson_vector *) F_PT(s,quark ), nrm );

  }


  /** FFT the convolution back to real space **/
  restrict_fourier_site(quark, sizeof(wilson_vector), 1); 

  /*** The 1 in the restrict_fourier routine corresponds to the 
       phase factor of exp( -i p x )
   ***/

  tt += dclock() ; 

  IF_VERBOSE_ON(1)
    printf("Time in sink_smear_light_pos_quark = %e sec\n",tt) ;

}




/*      >>>>>  spin_wilson_vector data structures <<<<<
 *   Smear the light quark propagator at the sink
 *
 *   quark_smear(x) = \sum_{r) quark(x+r) * smear(r)
 *
 *                  = \sum_{p) exp( -i p x)   quark(p) * smear(-p)
 *
 *   This function assumes that the quark propagator is in position space
 *   and that the smearing function is already in momentum space.
 *
 *
 *
 *
 *   Subroutine arguments 
 *
 *       quark        :: spin_wilson_vector pointer to the site structure, 
 *                     the quark propagator for a specific
 *                     source colour and spin in POSITION space.
 *                     The smeared quark propagator is returned in this 
 *                     data structure.
 *
 *       fft_work_a  :: site structure spin_wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *       fft_work_b  :: site structure spin_wilson_vector pointer 
 *                      [ work space for the FFT]
 *
 *      smear_func_mom :: site structure wilson_vector pointer 
 *                         smearing function in momentum space                
 *
 *   Any volume factors are multiplied into the correlators later in the code.
 *
 */


void sink_smear_light_pos_quark_large(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func_mom )
{
  register site *s;
  int i ; 
  double tt ;
  complex nrm ; 

  nrm.real = 1.0 / (nx*ny*nz)   ; nrm.imag = 0.0 ; 


  tt = -1.0*dclock() ; 

  /*** fourier transform the quark propagator 
    quark( p ) = sum_x exp( i p .x ) quark(x)

    This corresponds to the negative sign in restrict_fourier.
   *****/

  restrict_fourier_site(quark, sizeof(spin_wilson_vector), -1); 

  /** multiply the quark propagator by the smearing function in momentum space ****/
  FORALLSITES(i,s)
  {
    c_scalar_self_mult_spin_wvec((spin_wilson_vector *) F_PT(s,quark ),
			     *((complex *) F_PT(s, smear_func_mom)) );

    c_scalar_self_mult_spin_wvec((spin_wilson_vector *) F_PT(s,quark ), nrm ) ; 

  }


  /** FFT the convolution back to real space **/
  restrict_fourier_site(quark, sizeof(spin_wilson_vector), 1); 

  /*** The 1 in the restrict_fourier routine corresponds to the 
       phase factor of exp( -i p x )
   ***/

  tt += dclock() ; 

  IF_VERBOSE_ON(1)
    printf("Time in sink_smear_light_pos_quark = %e sec\n",tt) ;

}








