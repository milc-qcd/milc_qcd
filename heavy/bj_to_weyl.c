/************* bj_to_weyl.c   **************************/
/*
 * Multiply a Wilson vector (thought of as a column vector) by the  matrix V
 * acting on the left, where V converts between Weyl (w) and B&D (b)
 * conventions for the gamma matrices: gammab[mu] = Vadj gammaw[mu] V  . Note
 * that there is an extra minus sign for converting gamma5: gamma5b = - Vadj
 * gamma5w V usage:  bj_to_weyl( src, dest) wilson_vector *src,*dest; 
 *
 *
 * WARNING!!!  because the action convention has 1 + gamma_mu (as opposed to 1 -
 * gamma_mu) in the positive direction, it is the LOWER 2 components in B&D
 * conventions which hop forward!!!! 
 *
 *
Weyl conventions:
 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1


B&D conventions:
 gamma(XUP) 
 	    0  0  0 -i
            0  0 -i  0
            0  i  0  0
            i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0 -i  0
            0  0  0  i
            i  0  0  0
            0 -i  0  0

 gamma(TUP)
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1

 gamma(FIVE) 
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

Conversion matrix V
 	    0 -i  0  i
            i  0 -i  0    * (1/sqrt(2.))
            0 -i  0 -i    
            i  0  i  0


 */

#include "w_heavy_includes.h"


void bj_to_weyl(wilson_vector * src, wilson_vector * dest)
{
  register int i;		/* color */
  complex z1, z2 ;
  Real sqrt2inv;

  sqrt2inv = (Real) (1./ (sqrt(2.)));
  for (i = 0; i < 3; i++)
  {
    CSUB(src->d[3].c[i], src->d[1].c[i], z1);
    TIMESPLUSI(z1, z2);
    CMULREAL(z2, sqrt2inv, dest->d[0].c[i]);

    CSUB(src->d[0].c[i], src->d[2].c[i], z1);
    TIMESPLUSI(z1, z2);
    CMULREAL(z2, sqrt2inv, dest->d[1].c[i]);

    CADD(src->d[3].c[i], src->d[1].c[i], z1);
    TIMESMINUSI(z1, z2);
    CMULREAL(z2, sqrt2inv, dest->d[2].c[i]);

    CADD(src->d[0].c[i], src->d[2].c[i], z1);
    TIMESPLUSI(z1, z2);
    CMULREAL(z2, sqrt2inv, dest->d[3].c[i]);
  }
}
