/*
 *  This file contains a number of laguerre polnomials
 *
 *
 *
 */

#include<stdio.h>
#include<math.h>

/*
 * 2S  polynomial part of the hydrogen smearing functions
 */

Real two_s_hydrogen_poly(Real x)
{
  Real t0 ;

  t0 = 2.0*x-4.0;

  return t0 ;
}


/*
 * 3S polynomial part of the hydrogen smearing functions
 */

Real three_s_hydrogen_poly(Real x)
{
  Real t0 ;

  t0 = -18.0+(18.0-3.0*x)*x;

  return t0 ;
}







/*
 * 4S polynomial part of the hydrogen smearing functions
 */

Real four_s_hydrogen_poly(Real x)
{
  Real t0 ;

  t0 = -96.0+(144.0+(-48.0+4.0*x)*x)*x;

  return t0 ;
}








/*
 * 5S polynomial part of the hydrogen smearing functions
 */

Real five_s_hydrogen_poly(Real x)
{
  Real t0 ;

  t0 = -600.0+(1200.0+(-600.0+(100.0-5.0*x)*x)*x)*x;
  
  return t0 ;
}




