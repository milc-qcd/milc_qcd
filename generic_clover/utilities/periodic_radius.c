/*
 *  Calculate the radius on a periodic lattice
 *
 */

#include <stdio.h>
#include <math.h>
#include<string.h>



Real periodic_radius(int x, int y, int z, int nx)
{
  Real rad ;

  if( x >= nx/2 )  x -= nx ;
  if( y >= nx/2 )  y -= nx ;
  if( z >= nx/2 )  z -= nx ;

  rad = (Real) sqrt( (double)  x*x +  y*y +  z*z ) ;

  return rad ;

}



/**
 **  Allow for a squeeze in the x direction
 **
 **
 **/

Real periodic_radius_x_squeeze(int x, int y, int z, int nx, Real squeeze)
{
  Real rad ;

  if( x >= nx/2 )  x -= nx ;
  if( y >= nx/2 )  y -= nx ;
  if( z >= nx/2 )  z -= nx ;

  rad = (Real) sqrt( (double)  x*x*squeeze*squeeze  +  
		     y*y +  z*z ) ;

  return rad ;
}
