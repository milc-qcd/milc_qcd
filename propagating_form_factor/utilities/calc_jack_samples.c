/* 
 *   A collection of routines to calculate full and jackknife samples
 *   of the two and three point functions.
 *
 *
 *
 */


#include "read_hl_form.h"



/*   =====> two point   <=====
 *  Calculate the jackknife of full samples of the 
 *  two point correlators.
 *
 *     jomit :: sample to omit in the jackknife analysis
 *              [[ for full sample jomit = nosample)
 */

complex jack_sample_twopt(
  int jomit, 
  int nosample, 
  complex *two_corr[MAX_NO_FILE ] ,
  int t)
{
  Real norm ; 
  int i ; 
  complex ans = {0.0, 0.0} ; 

  assert( jomit >= 0 && jomit <= nosample ); 

  if( jomit != nosample )
    norm = nosample - 1 ; 
  else
    norm = nosample ; 


  for(i = 0 ; i < nosample ; ++i)
  {

    if( i != jomit )
      {
	ans.real += two_corr[i][t].real ; 
	ans.imag += two_corr[i][t].imag ; 
      }
  }

  ans.real /= norm;
  ans.imag /= norm;

  return ans ; 
}





/*    =======> three-point <==========
 *  Calculate the jackknife of full samples of the
 *  three point correlators.
 *
 *     jomit :: sample to omit in the jackknife analysis
 *              [[ for full sample jomit = nosample)
 */

complex jack_sample_threept(
  int jomit, 
  int nosample, 
  complex *three_corr[MAX_NO_FILE ] ,
  int t)
{
  Real norm ; 
  int i ; 
  complex ans = {0.0,0.0} ; 

  assert( jomit >= 0 && jomit <= nosample ); 

  if( jomit != nosample )
    norm = nosample - 1 ; 
  else
    norm = nosample ; 

  assert( norm != 0.0 ) ;


  for(i = 0 ; i < nosample ; ++i)
  {
    
    if( i != jomit )
      {
	ans.real += three_corr[i][t].real ; 
	ans.imag += three_corr[i][t].imag ; 
      }
  }

  ans.real /= norm;
  ans.imag /= norm;

  return ans ; 
}



