/*****************  grow4wvecs.c  (in su3.a) ****************************
*									*
*  If sum=0,								*
*  Grow and add four wilson_vectors 					*
*  If sum=1,								*
*  Grow and sum four wilson_vectors to another wilson_vector		*
* void grow_four_wvecs(a,b1,b2,b3,b4,sign,sum)				*
* wilson_vector *a; half_wilson_vector *b1,*b2,*b3,*b4;			*
* int sign,sum;								*
* A  <-  B1 + B2 + B3 + B4   or						*
* A  <-  A + B1 + B2 + B3 + B4						*
* B1 is expanded using gamma_x, B2 using gamma_y, etc. 			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"
/* grow and sum four wilson_vectors */

#ifndef FAST

void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum ){
    if(sum==0)wp_grow( b1,a,XUP,sign);
    else wp_grow_add( b1,a,XUP,sign);
    wp_grow_add( b2,a,YUP,sign);
    wp_grow_add( b3,a,ZUP,sign);
    wp_grow_add( b4,a,TUP,sign);
}

#else  /* "FAST" code has wp_grow_add inlined */
/* For the RS6000 */

void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum ){
  int i;
  if(sum==0)
    {
      /* wp_grow( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESMINUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESMINUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESPLUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESPLUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
    }
  else
    {
      /*wp_grow_add( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TMI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TMI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TPI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TPI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
    }
  
  /* wp_grow_add( b2,a,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUM( a->d[2].c[i], b2->h[1].c[i]);
	CSUB( a->d[3].c[i], b2->h[0].c[i], a->d[3].c[i] );
      }
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUB( a->d[2].c[i], b2->h[1].c[i], a->d[2].c[i] );
	CSUM( a->d[3].c[i], b2->h[0].c[i]);
      }
    }
  
  /* wp_grow_add( b3,a,ZUP,sign); */
  
  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TMI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TPI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  else
    {
      /* case ZDOWN:*/
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TPI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TMI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  
  /* wp_grow_add( b4,a,TUP,sign); */
  
  if(sign==PLUS)
    {
      /* case TUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b4->h[0].c[i]);
	CSUM( a->d[1].c[i], b4->h[1].c[i]);
	CSUM( a->d[2].c[i], b4->h[0].c[i]);
	CSUM( a->d[3].c[i], b4->h[1].c[i]);
      }
    }
  else
    {
      /* case TDOWN: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b4->h[0].c[i]);
	CSUM( a->d[1].c[i], b4->h[1].c[i]);
	CSUB( a->d[2].c[i], b4->h[0].c[i], a->d[2].c[i] );
	CSUB( a->d[3].c[i], b4->h[1].c[i], a->d[3].c[i] );
      }
    }
}

#endif /* "#ifndef FAST"*/
