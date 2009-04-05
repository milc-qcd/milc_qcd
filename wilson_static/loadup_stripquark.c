/******************** loadup_stripquark.c ******************************/
/* MIMD version 7 */
/*
 *  Routines to use the spin-colour source component
 *  of the quark propagator, in a wilson_vector data
 *  structure to build up the spinless object required in the 
 *  static calculation.
 */

#include "w_static_includes.h"

/*  trace( V (1 +/- gamma_4)) S_L )
 *
 * Conversion matrix V
 *           0 -i  0 -i
 *           i  0  i  0    * (1/sqrt(2.))
 *           0  i  0 -i
 *          -i  0  i  0
 *
 * The quark propagators use the Bjorken and Drell representation
 * at the source and the Weyl representation at the sink. The matrix V 
 * transforms the source represntation to Weyl.
 *
 *  This function assumes that the strip_quark array has been set
 *  to zero, before this routine has been called the first time.
 *
 *  Because of the nonstandard convention for the light quark
 *  action used by the MILC code.
 *
 *   strip_quark = trace( quark (1 - g4 ) )      t > 0
 *   strip_quark = trace( quark (1 + g4 ) )      t < 0
 *
 *
 */

void buildup_strip(field_offset quark, int colour, int spin)
{
  wilson_vector proj_quark ;
  wilson_vector g4_quark ;
  register site *s;
  int i ;
  int ic;

  complex mult_A , mult_B ;
  int spin_A = -1 , spin_B = -1;
  Real root_two= (Real) sqrt( (double) 2.0 ) ;
  complex z ;
  /*****------------------------------**********************/


  switch (spin) 
  {
  case 0 :
    spin_A = 1 ;
    mult_A = cmplx( 0.0 , -1.0/root_two ); 

    spin_B = 3 ;
    mult_B = cmplx( 0.0 , -1.0/root_two ); 

    break ;
  case 1 :
    spin_A = 0 ;
    mult_A = cmplx( 0.0 , 1.0/root_two ); 

    spin_B = 2 ;
    mult_B = cmplx( 0.0 , 1.0/root_two ); 

    break ;
  case 2 :

    spin_A = 1 ;
    mult_A = cmplx( 0.0 , 1.0/root_two ); 

    spin_B = 3 ;
    mult_B = cmplx( 0.0 , -1.0/root_two ); 

    break ;
  case 3 :

    spin_A = 0 ;
    mult_A = cmplx( 0.0 ,-1.0/root_two ); 

    spin_B = 2 ;
    mult_B = cmplx( 0.0 ,  1.0/root_two ); 

    break ;
  }

    FORALLSITES(i,s)
    {
      mult_by_gamma((wilson_vector *)F_PT(s,quark) , &(g4_quark), TUP );

      if( s->t < nt/2 ) 
	sub_wilson_vector((wilson_vector *)F_PT(s,quark), &(g4_quark)  ,&(proj_quark) );
      else
	add_wilson_vector((wilson_vector *)F_PT(s,quark), &(g4_quark)  ,&(proj_quark) );


       /** convert the source gamma matrix represenaion back to Weyl ****/
       for(ic=0 ; ic < 3 ;++ic)
       {

	CMUL(proj_quark.d[spin_A].c[ic] , mult_A , z) ;
        CSUM(s->strip_quark.e[ic][colour] , z ) ; 


	CMUL(proj_quark.d[spin_B].c[ic] , mult_B , z) ;
        CSUM(s->strip_quark.e[ic][colour] , z ) ; 

      }

    } /** end the loop over lattice sites ****/

}



