/**************************** create_smear_source.c ****************************/
/* MIMD version 7 */
/*
 *  This file contains a number of routines to calculate
 *  the smearing wave functions.
 *
 *
 */

#include "w_static_includes.h"

/*** functions required on this file ********/
Real periodic_radius(int x, int y, int z, int nx);

/*
 *  Create the local smearing operator.
 *
 *   src(x) = 0   if   x != 0
 *   src(x) = 1   if   x =  0
 *
 */


void create_local_oper(int which_smear, char filename[])
{
  int i ;
  register site *s;


    FORALLSITES(i,s)
    {
      if( s->x ==0 && s->y ==0 && s->z ==0 )
      {
	s->smear_func[which_smear].real = 1.0 ; 
	s->smear_func[which_smear].imag = 0.0 ; 
      }
      else
      {
	s->smear_func[which_smear].real = 0.0 ; 
	s->smear_func[which_smear].imag = 0.0 ; 
      }

    }  /** end loop over the lattice sites ****/




  sprintf(filename,"local_smear");

}




/*
 *  Create the wall smearing operator.
 *
 *   src(x) = 1   all x 
 *
 */


void create_wall_oper(int which_smear, char filename[])
{
  int i ;
  register site *s ;

    FORALLSITES(i,s)
    {
      s->smear_func[which_smear].real = 1.0 ; 
      s->smear_func[which_smear].imag = 0.0 ; 

    }  /** end loop over the lattice sites ****/


  sprintf(filename,"wall_smear");

}




/*
 *  Create an exponential smearing function
 *
 *   src(x) = exp( -alpha * x ) 
 *
 */


void create_expon_oper(int which_smear,   Real decay , char filename[])
{
  int i ;
  register site *s;
  Real rad ;


  /*** create the source ********/
  
  FORALLSITES(i,s)
  {
    rad =  periodic_radius(s->x, s->y, s->z, nx) ;
    
    s->smear_func[which_smear].real = (Real) exp( (double)  -decay*rad );
    s->smear_func[which_smear].imag = 0.0 ; 


  }  /*** end the loop over lattice sites ****/


  sprintf(filename,"exp_smear_%g",decay);

}


/*  2S model smearing function
 *
 *  Create an exponential times x smearing function
 *
 *   src(x) = exp( -alpha * x ) * (1 - A*x)
 *
 */


void create_2S_polyexpon_oper(int which_smear, Real decay, Real A, char filename[])
{
  Real rad ;
  int i ;
  register site *s;



  /*** create the source ********/
  
  FORALLSITES(i,s)
  {
    rad =  periodic_radius(s->x, s->y, s->z, nx) ;
    
    s->smear_func[which_smear].real = (Real) (1.0 - A*rad)*exp( (double)  -decay*rad );
    s->smear_func[which_smear].imag = 0.0 ; 

  }  /*** end the loop over the lattice ****/


  sprintf(filename,"polyexp_smear_%g_x%g",decay,A);


}



/*  3S model smearing function
 *
 *  Create an exponential times a polynomial smearing function
 *
 *   src(x) = exp( -alpha * x ) * (1 - B*x - D*x**2)
 *
 */


void create_3S_polyexpon_oper(int which_smear, Real decay, Real B, Real D,
char filename[])
{
  Real rad ;
  int i ;
  register site *s;



  /*** create the source ********/
  
  FORALLSITES(i,s)
  {
    
    rad =  periodic_radius(s->x, s->y, s->z, nx) ;
    
    s->smear_func[which_smear].real = 
      (Real) (1.0 - B*rad -D*rad*rad )*exp( (double)  -decay*rad );
    s->smear_func[which_smear].imag = 0.0 ; 

  }


  sprintf(filename,"polyexp_smear_%g_x%g_xsq%g",decay,B,D);

}



/*
 *  Calculate the radius on a periodic lattice
 *
 */

Real periodic_radius(int x, int y, int z, int dim)
{
  Real rad ;

  if( x >= dim/2 )  x -= dim ;
  if( y >= dim/2 )  y -= dim ;
  if( z >= dim/2 )  z -= dim ;

  rad = (Real) sqrt( (double)  x*x +  y*y +  z*z ) ;

  return rad ;

}
