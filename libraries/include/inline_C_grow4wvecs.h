/*****************  grow4wvecs.c  (in su3.a) ****************************
*									*
*  If sum=0,								*
*  Grow and add four wilson_vectors 					*
*  If sum=1,								*
*  Grow and sum four wilson_vectors to another wilson_vector		*
*  void grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,  *
*        half_wilson_vector *b2, half_wilson_vector *b3,                *
*	 half_wilson_vector *b4, int sign, int sum )                    *
* A  <-  B1 + B2 + B3 + B4   or						*
* A  <-  A + B1 + B2 + B3 + B4						*
* B1 is expanded using gamma_x, B2 using gamma_y, etc. 			*
*/
/* grow and sum four wilson_vectors */

#define _inline_C_grow_add_four_wvecs( a, b1, b2, b3, b4, sign, sum ) {\
  int _i;\
  if((sum)==0)\
    {\
      /* wp_grow( b1,a,XUP,sign); */\
      \
      /* case XUP: */\
      if((sign)==PLUS)\
	{\
	  for(_i=0;_i<3;_i++){\
	    (a)->d[0].c[_i]      = (b1)->h[0].c[_i];\
	    (a)->d[1].c[_i]      = (b1)->h[1].c[_i];\
	    TIMESMINUSI( (b1)->h[0].c[_i], (a)->d[3].c[_i]);\
	    TIMESMINUSI( (b1)->h[1].c[_i], (a)->d[2].c[_i]);\
	  }\
	}\
      else\
	{\
	  /* case XDOWN: */\
	  for(_i=0;_i<3;_i++){\
	    (a)->d[0].c[_i]      = (b1)->h[0].c[_i];\
	    (a)->d[1].c[_i]      = (b1)->h[1].c[_i];\
	    TIMESPLUSI( (b1)->h[0].c[_i], (a)->d[3].c[_i]);\
	    TIMESPLUSI( (b1)->h[1].c[_i], (a)->d[2].c[_i]);\
	  }\
	}\
    }\
  else\
    {\
      /*wp_grow_add( b1,a,XUP,sign); */\
      \
      /* case XUP: */\
      if((sign)==PLUS)\
	{\
	  for(_i=0;_i<3;_i++){\
	    CSUM( (a)->d[0].c[_i], (b1)->h[0].c[_i]);\
	    CSUM( (a)->d[1].c[_i], (b1)->h[1].c[_i]);\
	    CSUM_TMI( (a)->d[2].c[_i], (b1)->h[1].c[_i] );\
	    CSUM_TMI( (a)->d[3].c[_i], (b1)->h[0].c[_i] );\
	  }\
	}\
      else\
	{\
	  /* case XDOWN: */\
	  for(_i=0;_i<3;_i++){\
	    CSUM( (a)->d[0].c[_i], (b1)->h[0].c[_i]);\
	    CSUM( (a)->d[1].c[_i], (b1)->h[1].c[_i]);\
	    CSUM_TPI( (a)->d[2].c[_i], (b1)->h[1].c[_i] );\
	    CSUM_TPI( (a)->d[3].c[_i], (b1)->h[0].c[_i] );\
	  }\
	}\
    }\
  \
  /* wp_grow_add( b2,a,YUP,sign); */\
  \
  if((sign)==PLUS)\
    {\
      /* case YUP: */\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b2)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b2)->h[1].c[_i]);\
	CSUM( (a)->d[2].c[_i], (b2)->h[1].c[_i]);\
	CSUB( (a)->d[3].c[_i], (b2)->h[0].c[_i], (a)->d[3].c[_i] );\
      }\
    }\
  else\
    {\
      /* case YDOWN: */\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b2)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b2)->h[1].c[_i]);\
	CSUB( (a)->d[2].c[_i], (b2)->h[1].c[_i], (a)->d[2].c[_i] );\
	CSUM( (a)->d[3].c[_i], (b2)->h[0].c[_i]);\
      }\
    }\
  \
  /* wp_grow_add( b3,a,ZUP,sign); */\
  \
  if((sign)==PLUS)\
    {\
      /* case ZUP: */\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b3)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b3)->h[1].c[_i]);\
	CSUM_TMI( (a)->d[2].c[_i], (b3)->h[0].c[_i] );\
	CSUM_TPI( (a)->d[3].c[_i], (b3)->h[1].c[_i] );\
      }\
    }\
  else\
    {\
      /* case ZDOWN:*/\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b3)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b3)->h[1].c[_i]);\
	CSUM_TPI( (a)->d[2].c[_i], (b3)->h[0].c[_i] );\
	CSUM_TMI( (a)->d[3].c[_i], (b3)->h[1].c[_i] );\
      }\
    }\
  \
  /* wp_grow_add( b4,a,TUP,sign); */\
  \
  if((sign)==PLUS)\
    {\
      /* case TUP: */\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b4)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b4)->h[1].c[_i]);\
	CSUM( (a)->d[2].c[_i], (b4)->h[0].c[_i]);\
	CSUM( (a)->d[3].c[_i], (b4)->h[1].c[_i]);\
      }\
    }\
  else\
    {\
      /* case TDOWN: */\
      for(_i=0;_i<3;_i++){\
	CSUM( (a)->d[0].c[_i], (b4)->h[0].c[_i]);\
	CSUM( (a)->d[1].c[_i], (b4)->h[1].c[_i]);\
	CSUB( (a)->d[2].c[_i], (b4)->h[0].c[_i], (a)->d[2].c[_i] );\
	CSUB( (a)->d[3].c[_i], (b4)->h[1].c[_i], (a)->d[3].c[_i] );\
      }\
    }\
}
