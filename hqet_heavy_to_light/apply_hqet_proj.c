/**************************** apply_hqet_proj.c ****************************/
/* MIMD version 6 */
/*
   I have HACKED the code so that I get the 
   correct convention for the upper and 
   lowere components.



   Apply the HQET projection operator
   to a quark propagator

  t <= tf     quark = P_[-](v) * quark  (backwards moving) 
  t > tf      quark = P_[+](v) * quark  (forwards moving) 


  P_[+](v) = 1/2( 1 + v_0 * \gamma_0 - i *  v_k * \gamma_k )

  P_[-](v) = 1/2( 1 - v_0 * \gamma_0 + i *  v_k * \gamma_k )


  Program notes


*  Multiply a Wilson vector by a  complex scalar and add to another vector
* dest  <-  src1 + s*src2


 */

/* MIMD version 6 */

/* Initialize a source for the inverter */

#include "hqet_light_includes.h"

/*
 *
 *  Function arguments 
 *    on input
 *       ans   ::  quark propagator in a wilson vector pointer
 *                 to the site structure
 *        v_pt ::  which velocity to use
 *
 *    on output
 *       ans   :: the quark propagator multiplied by the projection
 *                operator.
 *
 *
 */

void apply_hqet_proj(field_offset ans, int v_pt)
{
  register int i;
  register site *s; 
  complex z ; 
  wilson_vector proj ; 
  wilson_vector tmp ; 
  Real sgn = 1 ; 


if( this_node == 0 ) printf("MILC quark --> g5 quark g5 convention corrected\n") ; 


  FORALLSITES(i,s) 
  {
    if( s->t < tf  )
    {
      sgn = 1.0 ; 
    }
    else
    {
     sgn = -1.0 ; 
    }

    clear_wvec(&proj);

    /**** unit matrix contribution ****/
    z = cmplx( 0.5 , 0.0 ) ;
    c_scalar_mult_add_wvec(&proj, (wilson_vector *)F_PT(s,ans),&z , &proj) ; 

    /******   gamma_time *********/
    mult_by_gamma((wilson_vector *)F_PT(s,ans), &tmp, TUP );
    z = cmplx(sgn*0.5*velocity[ v_pt ][ TUP] , 0.0 ) ; 
    c_scalar_mult_add_wvec(&proj, &tmp , &z , &proj) ; 


    /******   gamma_X *********/
    mult_by_gamma((wilson_vector *)F_PT(s,ans), &tmp, XUP );
    z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ XUP] ) ;
    c_scalar_mult_add_wvec(&proj, &tmp , &z , &proj) ; 
		    
    /******   gamma_Y *********/
    mult_by_gamma((wilson_vector *)F_PT(s,ans), &tmp, YUP );
    z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ YUP ] ) ; 
    c_scalar_mult_add_wvec(&proj, &tmp , &z , &proj) ; 
		    
    /******   gamma_Z *********/
    mult_by_gamma((wilson_vector *)F_PT(s,ans), &tmp, ZUP );
    z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ ZUP ] ) ;
    c_scalar_mult_add_wvec(&proj, &tmp , &z , &proj) ; 


    /*** copy the projection operator into the orginal source ***/
    *((wilson_vector *)F_PT(s,ans)) = proj ; 
    


  } /*** end the loop over sites on the lattice ****/



}  /*** end of the function *****/






