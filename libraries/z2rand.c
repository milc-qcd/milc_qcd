/*********** z2rand.c (in su3.a) ****************************************
*									*
*  int z2_ran_no( double_prn *prn_pt )   				*
*  Z(2)/sqrt(2) distributed random number 				*
*  This requires a random number generator named "myrand()", returning	*
*  a Real uniformly distributed between zero and one. The argument of	*
*  this routine is a pointer to be passed to myrand(). 			*
*/

#include "../include/config.h"
#include <math.h>
#include "../include/su3.h"
#include "../include/random.h"

Real z2_rand_no( double_prn *prn_pt ){

  if( myrand(prn_pt) > 0.5 ) return 1./sqrt(2.);
  else return -1./sqrt(2.);
}
