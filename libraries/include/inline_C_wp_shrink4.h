/* Included here to allow for inlining */
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

#define _inline_C_wp_shrink_4dir( a, b1, b2, b3, b4, sign ) {\
  register int _i; /*color*/\
  \
/*    wp_shrink( a,b1,XUP,sign); */\
\
  if((sign)==PLUS)\
    {\
    /* case XUP: */\
      for(_i=0;_i<3;_i++){\
	(b1)->h[0].c[_i].real = (a)->d[0].c[_i].real - (a)->d[3].c[_i].imag;\
	(b1)->h[0].c[_i].imag = (a)->d[0].c[_i].imag + (a)->d[3].c[_i].real;\
	(b1)->h[1].c[_i].real = (a)->d[1].c[_i].real - (a)->d[2].c[_i].imag;\
	(b1)->h[1].c[_i].imag = (a)->d[1].c[_i].imag + (a)->d[2].c[_i].real;\
      }\
    }\
  else\
    {\
      /* case XDOWN: */\
      for(_i=0;_i<3;_i++){\
	(b1)->h[0].c[_i].real = (a)->d[0].c[_i].real + (a)->d[3].c[_i].imag;\
	(b1)->h[0].c[_i].imag = (a)->d[0].c[_i].imag - (a)->d[3].c[_i].real;\
	(b1)->h[1].c[_i].real = (a)->d[1].c[_i].real + (a)->d[2].c[_i].imag;\
	(b1)->h[1].c[_i].imag = (a)->d[1].c[_i].imag - (a)->d[2].c[_i].real;\
      }\
    }\
  \
  \
  /*    wp_shrink( a,b2,YUP,sign); */\
  \
  if((sign)==PLUS)\
    {\
      /* case YUP: */\
      for(_i=0;_i<3;_i++){\
	(b2)->h[0].c[_i].real = (a)->d[0].c[_i].real - (a)->d[3].c[_i].real;\
	(b2)->h[0].c[_i].imag = (a)->d[0].c[_i].imag - (a)->d[3].c[_i].imag;\
	(b2)->h[1].c[_i].real = (a)->d[1].c[_i].real + (a)->d[2].c[_i].real;\
	(b2)->h[1].c[_i].imag = (a)->d[1].c[_i].imag + (a)->d[2].c[_i].imag;\
      }\
      \
    }\
  else\
    {\
      /* case YDOWN: */\
      for(_i=0;_i<3;_i++){\
	(b2)->h[0].c[_i].real = (a)->d[0].c[_i].real + (a)->d[3].c[_i].real;\
	(b2)->h[0].c[_i].imag = (a)->d[0].c[_i].imag + (a)->d[3].c[_i].imag;\
	(b2)->h[1].c[_i].real = (a)->d[1].c[_i].real - (a)->d[2].c[_i].real;\
	(b2)->h[1].c[_i].imag = (a)->d[1].c[_i].imag - (a)->d[2].c[_i].imag;\
      }\
    }\
  \
  /*    wp_shrink( a,b3,ZUP,sign); */\
\
  if((sign)==PLUS)\
    {\
      /* case ZUP: */\
      for(_i=0;_i<3;_i++){\
	(b3)->h[0].c[_i].real = (a)->d[0].c[_i].real - (a)->d[2].c[_i].imag;\
	(b3)->h[0].c[_i].imag = (a)->d[0].c[_i].imag + (a)->d[2].c[_i].real;\
	(b3)->h[1].c[_i].real = (a)->d[1].c[_i].real + (a)->d[3].c[_i].imag;\
	(b3)->h[1].c[_i].imag = (a)->d[1].c[_i].imag - (a)->d[3].c[_i].real;\
      }\
    }\
  else\
    {\
      /* case ZDOWN: */\
      for(_i=0;_i<3;_i++){\
	(b3)->h[0].c[_i].real = (a)->d[0].c[_i].real + (a)->d[2].c[_i].imag;\
	(b3)->h[0].c[_i].imag = (a)->d[0].c[_i].imag - (a)->d[2].c[_i].real;\
	(b3)->h[1].c[_i].real = (a)->d[1].c[_i].real - (a)->d[3].c[_i].imag;\
	(b3)->h[1].c[_i].imag = (a)->d[1].c[_i].imag + (a)->d[3].c[_i].real;\
      }\
      \
    }\
\
/*    wp_shrink( a,b4,TUP,sign); */\
\
  if((sign)==PLUS)\
    {\
      /* case TUP: */\
      for(_i=0;_i<3;_i++){\
	(b4)->h[0].c[_i].real = (a)->d[0].c[_i].real + (a)->d[2].c[_i].real;\
	(b4)->h[0].c[_i].imag = (a)->d[0].c[_i].imag + (a)->d[2].c[_i].imag;\
	(b4)->h[1].c[_i].real = (a)->d[1].c[_i].real + (a)->d[3].c[_i].real;\
	(b4)->h[1].c[_i].imag = (a)->d[1].c[_i].imag + (a)->d[3].c[_i].imag;\
      }\
    }\
  else\
    {\
      /* case TDOWN: */\
      for(_i=0;_i<3;_i++){\
	(b4)->h[0].c[_i].real = (a)->d[0].c[_i].real - (a)->d[2].c[_i].real;\
	(b4)->h[0].c[_i].imag = (a)->d[0].c[_i].imag - (a)->d[2].c[_i].imag;\
	(b4)->h[1].c[_i].real = (a)->d[1].c[_i].real - (a)->d[3].c[_i].real;\
	(b4)->h[1].c[_i].imag = (a)->d[1].c[_i].imag - (a)->d[3].c[_i].imag;\
      }\
    }\
}
