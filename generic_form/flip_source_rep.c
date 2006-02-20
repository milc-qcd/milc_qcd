/**************** flip_source_rep.c ********************************/  
/* MIMD version 7 */
/* Flip the gamma matrix represenatation at the
  source. This function is required to use 
  quark propagators from Claude's hopping variational code.

 ***/


#include "generic_form_includes.h"

/* function prototypes */
void define_change_rep(complex change_rep[4][4] ) ;
void gammafive_sandwich(field_offset prop) ;

/*************************************************************

  quark_prop = V  quark_prop

  where V is a gamma matrix, that operates on the source
  spin index.

   Function arguments
   -------------------

     On input 
        quark_prop ::  site structure pointer to the 
                       quark propagator in a spin_wilson_vector file.

     On output
        quark_prop ::  site structure pointer to the 
                       quark propagator in a spin_wilson_vector file.


***************************************************************/


void flip_source_re(field_offset quark_prop)
{
  register int i ;
  register site *s; 
  static int flag_call = 0 ;
  
  complex change_rep[4][4] ;

  complex z ; 
  spin_wilson_vector quark_flip;          /* temporary storage for quark */

  int si ,  cf , sf ;
  int k ; 


  /**** trick to print out some titles only on the first call ****/
  if( flag_call == 0 )
  {
    IF_MASTER
      printf("\nRotation of the source gamma representation has been set up\n");

    flag_call = 1 ;
  }


  define_change_rep( change_rep ); 

  FORALLSITES(i,s)
  {

	for(si=0;si<4;si++)
	  for(sf=0;sf<4;sf++)
	    for(cf=0;cf<3;cf++)
	    {
	      quark_flip.d[si].d[sf].c[cf].real = 0.0  ; 
	      quark_flip.d[si].d[sf].c[cf].imag = 0.0  ; 

	      for(k = 0 ; k < 4  ; ++k)
	      {
/**		CMUL(((spin_wilson_vector *)F_PT(s,quark_prop))->d[k].d[sf].c[cf] , change_rep[sf][k] ,z)  ;  *** still wrong ***/
		CMUL_J(((spin_wilson_vector *)F_PT(s,quark_prop))->d[k].d[sf].c[cf] , change_rep[si][k] ,z)  ; 
		CSUM(quark_flip.d[si].d[sf].c[cf]  ,z) ; 
	      }

	    }  /*** end the loop over colour and spin *****/

	*((spin_wilson_vector *)F_PT(s,quark_prop)) =   quark_flip ; 


      }  /**** end of the loop over lattice sites ******/
  

  /**** flip the action convention for forwards and backwards *****/
/**   gammafive_sandwich(quark_prop) ;   ***/



} /* end of the flip function  */




/*
 *  Define the matrix that flipe the representation of the 
 *  source gamma matrices.
 *


*
  Multiply a Wilson vector (thought of as a column vector)
   by the  matrix V acting on the left,
   where V converts between Weyl (w) and B&D (b)
   conventions for the gamma matrices:
   gammab[mu] = Vadj gammaw[mu] V  .
   Note that there is an extra minus sign
   for converting gamma5:
   gamma5b = - Vadj gamma5w V
     usage:  bj_to_weyl( src, dest)
        wilson_vector *src,*dest;


        WARNING!!!  because the action convention has 1 + gamma_mu
        (as opposed to 1 - gamma_mu) in the positive direction, it is
        the LOWER 2 components in B&D conventions which hop forward!!!!


Conversion matrix V
            0 -i  0  i
            i  0 -i  0    * (1/sqrt(2.))
            0 -i  0 -i
            i  0  i  0

 *
 */


void define_change_rep(complex change_rep[4][4] )
{
  int i , j ;
  Real fact = (Real) 1.0/sqrt(2.0) ; 


  for(i=0 ; i < 4 ; ++i)
    for(j=0 ; j < 4 ; ++j)
    {
      change_rep[i][j].real = 0.0 ; 
      change_rep[i][j].imag = 0.0 ; 
    }



  change_rep[0][1].imag = -fact  ; 
  change_rep[0][3].imag =  fact  ; 

  change_rep[1][0].imag =  fact  ; 
  change_rep[1][2].imag = -fact  ; 

  change_rep[2][1].imag = -fact  ; 
  change_rep[2][3].imag = -fact  ; 

  change_rep[3][0].imag =  fact  ; 
  change_rep[3][2].imag =  fact  ; 


#define DISPLAY_CHANGE_REP 1

#ifdef DISPLAY_CHANGE_REP

  IF_MASTER
  {

    printf("Here is the gamma matrix conversion matrix \n");

    for(i =0 ; i < 4 ; ++i)
    {
      printf("[ ");
      for(j =0 ; j < 4 ; ++j)
      {
	printf(" (%g,%g) ",change_rep[i][j].real , change_rep[i][j].imag ); 
      }
      printf(" ]\n");
    }

  }

#endif

}



/*************************************************************
  Flip the quark propagators into a standard representatiom

   prop -->  gamma_5  prop gamma_5


This transformation is requires because of Claude's comment


	WARNING!!!  because the action convention has 1 + gamma_mu
	(as opposed to 1 - gamma_mu) in the positive direction, it is
	the LOWER 2 components in B&D conventions which hop forward!!!!


Subroutine arguments
   prop  ::  site structure pouinter to spin_wilson_vector,
             which contains the quark propagator for a specific
             color.
  
************************************************************/

void gammafive_sandwich(field_offset prop)
{
  register int i;
  register site *s; 
  spin_wilson_vector localmat;       /* temporary storage for quark */
  spin_wilson_vector quark;       /* temporary storage for quark */
  int si,sf, cf;

  FORALLSITES(i,s)
  {
    mult_sw_by_gamma_l( (spin_wilson_vector *)F_PT(s,prop) , &localmat, GAMMAFIVE);    
    mult_sw_by_gamma_r( &localmat, &quark, GAMMAFIVE);    

    for(si=0;si<4;si++)
      for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++)
	{
	  ((spin_wilson_vector *)F_PT(s,prop))->d[si].d[sf].c[cf] = 
	   quark.d[si].d[sf].c[cf] ;
	    
	}




  }

}




