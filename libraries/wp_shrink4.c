/*****************  wp_shrink4.c  (in su3.a) ****************************
*									*
* Shrink a wilson vector in four directions, producing four		*
*  half_wilson_vectors.							*
* void wp_shrink_4dir(  wilson_vector *a,  half_wilson_vector *b1,	*
*       half_wilson_vector *b2, half_wilson_vector *b3,			*
*       half_wilson_vector *b4, int sign );				*
* B1 <- (1 +- gamma_x)A,, projection					*
*  argument "sign" is sign of gamma matrix.				*
*  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

#ifndef FAST   /* "FAST", or IBM RS6000 version inlines calls */

void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
        half_wilson_vector *b4, int sign ){
    wp_shrink( a,b1,XUP,sign);
    wp_shrink( a,b2,YUP,sign);
    wp_shrink( a,b3,ZUP,sign);
    wp_shrink( a,b4,TUP,sign);
}

#else   /* "FAST" code inlines calls */

void wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
        half_wilson_vector *b4, int sign ){
  register int i; /*color*/
  
/*    wp_shrink( a,b1,XUP,sign); */

  if(sign==PLUS)
    {
    /* case XUP: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].real;
      }
    }
  else
    {
      /* case XDOWN: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].real;
      }
    }
  
  
  /*    wp_shrink( a,b2,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].imag;
      }
      
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].imag;
      }
    }
  
  /*    wp_shrink( a,b3,ZUP,sign); */

  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real - a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag + a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real + a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag - a->d[3].c[i].real;
      }
    }
  else
    {
      /* case ZDOWN: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real + a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag - a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real - a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag + a->d[3].c[i].real;
      }
      
    }

/*    wp_shrink( a,b4,TUP,sign); */

  if(sign==PLUS)
    {
      /* case TUP: */
      for(i=0;i<3;i++){
	b4->h[0].c[i].real = a->d[0].c[i].real + a->d[2].c[i].real;
	b4->h[0].c[i].imag = a->d[0].c[i].imag + a->d[2].c[i].imag;
	b4->h[1].c[i].real = a->d[1].c[i].real + a->d[3].c[i].real;
	b4->h[1].c[i].imag = a->d[1].c[i].imag + a->d[3].c[i].imag;
      }
    }
  else
    {
      /* case TDOWN: */
      for(i=0;i<3;i++){
	b4->h[0].c[i].real = a->d[0].c[i].real - a->d[2].c[i].real;
	b4->h[0].c[i].imag = a->d[0].c[i].imag - a->d[2].c[i].imag;
	b4->h[1].c[i].real = a->d[1].c[i].real - a->d[3].c[i].real;
	b4->h[1].c[i].imag = a->d[1].c[i].imag - a->d[3].c[i].imag;
      }
    }
}

#endif /* "ifndef FAST */
