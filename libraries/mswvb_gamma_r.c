/************* mswvb_gamma_r.c **************************/
/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_swv_by_gamma_r( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir)
{
register int i; /*color*/
register int s1;	/* row  spin indices*/

  switch(dir){
    case XUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSI( src->d[s1].d[3].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESMINUSI( src->d[s1].d[2].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].d[1].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESPLUSI(  src->d[s1].d[0].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case YUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[2].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[1].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].d[0].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSI( src->d[s1].d[2].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSI(  src->d[s1].d[3].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].d[0].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSI( src->d[s1].d[1].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESPLUSONE( src->d[s1].d[2].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSONE( src->d[s1].d[0].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESPLUSONE( src->d[s1].d[1].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESPLUSONE(  src->d[s1].d[0].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[1].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESMINUSONE( src->d[s1].d[2].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
  }
}
