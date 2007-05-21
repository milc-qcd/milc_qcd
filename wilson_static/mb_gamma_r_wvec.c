/************* mb_gamma_r_wvec.c  (in su3.a) **************************/
/* MIMD version 7 */
/* 
  Multiply a Wilson vector (thought of as a row vector)
   by a gamma matrix acting on the right 
  usage:  mult_by_gamma_right( src, dest, dir)
	wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE

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
*/
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_by_gamma_right_vec( wilson_vector *src, wilson_vector *dest, int dir)
{
register int i; /*color*/

  switch(dir){
    case XUP:
	for(i=0;i<3;i++){
	    TIMESMINUSI( src->d[3].c[i], dest->d[0].c[i] );
	    TIMESMINUSI( src->d[2].c[i], dest->d[1].c[i] );
	    TIMESPLUSI(  src->d[1].c[i], dest->d[2].c[i] );
	    TIMESPLUSI(  src->d[0].c[i], dest->d[3].c[i] );
	}
	break;
    case YUP:
	for(i=0;i<3;i++){
	    TIMESMINUSONE( src->d[3].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE(  src->d[2].c[i], dest->d[1].c[i] );
	    TIMESPLUSONE(  src->d[1].c[i], dest->d[2].c[i] );
	    TIMESMINUSONE( src->d[0].c[i], dest->d[3].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<3;i++){
	    TIMESMINUSI( src->d[2].c[i], dest->d[0].c[i] );
	    TIMESPLUSI(  src->d[3].c[i], dest->d[1].c[i] );
	    TIMESPLUSI(  src->d[0].c[i], dest->d[2].c[i] );
	    TIMESMINUSI( src->d[1].c[i], dest->d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<3;i++){
	    TIMESPLUSONE( src->d[2].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE( src->d[3].c[i], dest->d[1].c[i] );
	    TIMESPLUSONE( src->d[0].c[i], dest->d[2].c[i] );
	    TIMESPLUSONE( src->d[1].c[i], dest->d[3].c[i] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<3;i++){
	    TIMESPLUSONE(  src->d[0].c[i], dest->d[0].c[i] );
	    TIMESPLUSONE(  src->d[1].c[i], dest->d[1].c[i] );
	    TIMESMINUSONE( src->d[2].c[i], dest->d[2].c[i] );
	    TIMESMINUSONE( src->d[3].c[i], dest->d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_RIGHT_VEC()\n");
  }
}

