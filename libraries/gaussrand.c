/*****************  gaussrand.c  (in su3.a) *****************************
*									*
*  Real gaussian_ran_no( double_prn *prn_pt )				*
*  Gaussian distributed random number					*
*  Probability distribution exp( -x*x ), so < x^2 > = 1/2		*
*  This requires a random number generator named "myrand()", returning	*
*  a Real uniformly distributed between zero and one. The argument of	*
*  this routine is a pointer to be passed to myrand(). 			*
*/

#include "../include/config.h"
#include <math.h>
#include "../include/su3.h"
#include "../include/random.h"

Real gaussian_rand_no( double_prn *prn_pt ){
static int iset=0;
static Real gset;
Real fac,r,v1,v2;

    if  (iset == 0) {
	do {
	    v1=2.0*myrand(prn_pt)-1.0;
	    v2=2.0*myrand(prn_pt)-1.0;
	    r=v1*v1+v2*v2;
	} while (r >= 1.0);
	fac=sqrt( -log((double)r)/(double)r);
	gset=v1*fac;
	iset=1;
	return v2*fac;
    } else {
	iset=0;
	return gset;
    }
}
