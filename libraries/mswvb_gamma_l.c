/************* mswvb_gamma_l.c **************************/
/* 
  Multiply a spin_wilson_vector by a gamma matrix
  acting on the row index
  (This is the first index, or equivalently, multiplication on the left)
  usage:  mb_gamma_l( src, dest, dir)
	spin_wilson_vector *src,*dest;
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
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_swv_by_gamma_l( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir)
{
register int c2,s2;	/* column indices, color and spin */

  switch(dir){
    case XUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSI(  src->d[3].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[2].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSI( src->d[1].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case YUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESMINUSONE( src->d[3].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[2].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[0].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case ZUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSI(  src->d[2].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESMINUSI( src->d[3].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[1].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case TUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSONE( src->d[2].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[3].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[0].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[1].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case GAMMAFIVE:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSONE(  src->d[0].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[2].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[3].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_LEFT()\n");
  }
}
