/*
 *   A set of routines to calculate jacknife errors.
 *
 *
 *
 * 
 *
 */

#include<math.h>
#include<stdio.h>
#include<string.h>



/*  A function to calculate the mean possibly leaving out
 *  an element.
 *
 *  data[0..nodata-1] X conventions.
 * 
 *  Single elimination jackknife.
 *
 *  if jomit = 0 then the standard average is taken.
 *
 */   
Real jackmean(Real *data,int nodata,int jomit)
{
  int i;
  Real mean = 0.0;

  for(i=0;i< nodata;++i)
    if( i != jomit ) mean += data[i] ;


  if( jomit != nodata )
    mean /= 1.0*(nodata-1);
  else
    mean /= 1.0*(nodata);


  return mean;

}


/*
 *  Calculate the error in a table of values using the
 *  jackknife method.
 *
 *  data[0..nodata-1] C conventions.
 *
 */

Real jackerror(Real mean,Real *data,int nodata)
{
  int i;
  Real error=0.0;


  for(i=0;i< nodata;++i)
    error += pow(mean - data[i],2.0);


  return (Real) sqrt( (double)  (nodata*1.0-1.0)*error/(1.0*nodata ));

}




