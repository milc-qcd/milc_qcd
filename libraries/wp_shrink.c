/************* wp_shrink.c  (in su3.a) **************************/
/* 
  Compute the "Wilson projection" of a Wilson fermion vector.
  (1 +- gamma_j) is a projection operator, and we want to isolate
  the components of the vector that it keeps.  In other words, keep
  the components of the vector along the eigenvectors of 1+-gamma_j
  with eigenvalue 2, and throw away those with eigenvalue 0.

  usage:  wp_shrink( wilson_vector *src, half_wilson_vector *dest,
	int dir, int sign )

	If dir is one of XUP,YUP,ZUP or TUP, take the projections
	along the eigenvectors with eigenvalue +1, which survive
	multiplication by (1+gamma[dir]).
	If dir is one of XDOWN,YDOWN,ZDOWN or TDOWN, take the projections
	along the eigenvectors with eigenvalue -1, which survive
	multiplication by (1-gamma[OPP_DIR(dir)]).
	If sign=MINUS, switch the roles of +1 and -1 (ie use -gamma_dir
	instead of gamma_dir )

  Here my eigenvectors are normalized to 2, so for XYZT directions
  I won't explicitely multiply by 2.  In other words, the matrix of
  eigenvectors is sqrt(2) times a unitary matrix, and in reexpanding
  the vector I will multiply by the adjoint of this matrix.

  For UP directions, hvec.h[0] and hvec.h[2] contain the projections
  along the first and second eigenvectors respectively.
  For DOWN directions, hvec.h[0] and hvec.h[2] contain the projections
  along the third and fourth eigenvectors respectively. This results
  in down directions differing from up directions only in the sign of
  the addition.

  Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to multiplication
   by 1+-gamma_dir

 gamma(XUP) 			eigenvectors	eigenvalue
 	    0  0  0  i		( 1, 0, 0,-i)	+1
            0  0  i  0		( 0, 1,-i, 0)	+1
            0 -i  0  0		( 0, 1, 0,+i)	-1
           -i  0  0  0		( 1, 0,+i, 0)	-1

 gamma(YUP)			eigenvectors	eigenvalue
 	    0  0  0 -1		( 1, 0, 0,-1)	+1
            0  0  1  0		( 0, 1, 1, 0)	+1
            0  1  0  0		( 1, 0, 0, 1)	-1
           -1  0  0  0		( 0, 1,-1, 0)	-1

 gamma(ZUP)			eigenvectors	eigenvalue
 	    0  0  i  0		( 1, 0,-i, 0)	+1
            0  0  0 -i		( 0, 1, 0,+i)	+1
           -i  0  0  0		( 1, 0,+i, 0)	-1
            0  i  0  0		( 0, 1, 0,-i)	-1

 gamma(TUP)			eigenvectors	eigenvalue
 	    0  0  1  0		( 1, 0, 1, 0)	+1
            0  0  0  1		( 0, 1, 0, 1)	+1
            1  0  0  0		( 1, 0,-1, 0)	-1
            0  1  0  0		( 0, 1, 0,-1)	-1

 gamma(FIVE) 			eigenvectors	eigenvalue
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void wp_shrink( wilson_vector *src, half_wilson_vector *dest,
        int dir, int sign ){
  register int i; /*color*/

  if(sign==MINUS)dir=OPP_DIR(dir);	/* two ways to get -gamma_dir ! */
  switch(dir){
    case XUP:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[3].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[3].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[2].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[2].c[i].real;
	}
	break;
    case XDOWN:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[3].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[3].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[2].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[2].c[i].real;
	}
	break;
    case YUP:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[3].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[3].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[2].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[2].c[i].imag;
	}
	break;
    case YDOWN:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[3].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[3].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[2].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[2].c[i].imag;
	}
	break;
    case ZUP:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[2].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[2].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[3].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[3].c[i].real;
	}
	break;
    case ZDOWN:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[2].c[i].imag;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[2].c[i].real;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[3].c[i].imag;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[3].c[i].real;
	}
	break;
    case TUP:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real + src->d[2].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag + src->d[2].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real + src->d[3].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag + src->d[3].c[i].imag;
	}
	break;
    case TDOWN:
	for(i=0;i<3;i++){
	    dest->h[0].c[i].real = src->d[0].c[i].real - src->d[2].c[i].real;
	    dest->h[0].c[i].imag = src->d[0].c[i].imag - src->d[2].c[i].imag;
	    dest->h[1].c[i].real = src->d[1].c[i].real - src->d[3].c[i].real;
	    dest->h[1].c[i].imag = src->d[1].c[i].imag - src->d[3].c[i].imag;
	}
	break;
    default:
	printf("BAD CALL TO WP_SHRINK()\n");
  }
}

