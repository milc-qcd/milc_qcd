/***************** wp_grow_pl_field_pf_l.c  **************/
/* 
  Expand the "Wilson projection" of a Wilson fermion vector.
  (1 +- gamma_j) is a projection operator, and we are given a
  half_wilson_vector which contains the two components of a Wilson
  vector projected out.  This routine reexpands it to a four component
  object.

  usage:  wp_grow(  half_wilson_vector *src, wilson_vector *dest,
        int dir, int sign );



  Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to multiplication
   by 1+-gamma_dir

 ahead    0: 1 1 0 0      8:
 ahead    1: 1 0 1 0      10:
 ahead    2: 1 0 0 1      12:
 ahead    3: 1 0 0 0      14:
 ahead    4: 1 0 0 -1      16:
 ahead    5: 1 0 -1 0      18:
 ahead    6: 1 -1 0 0      20:
 ahead    7: 0 1 1 0      22:
 ahead    8: 0 1 0 1      24:
 ahead    9: 0 1 0 0      26:
 ahead    10: 0 1 0 -1      28:
 ahead    11: 0 1 -1 0      30:
 ahead    12: 0 0 1 1      32:
 ahead    13: 0 0 1 0      34:
 ahead    14: 0 0 1 -1      36:
 ahead    15: 0 0 0 1      38:

*/
#include "arb_ov_includes.h"
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

#define SR05 0.707106781187
#define SR08 0.353553390593


#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif

#ifdef PREFETCH
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#endif


void wp_grow_pl_field_l(  half_wilson_vector * __restrict__ src, wilson_vector * __restrict__ dest,
        const int dir, const int sign , const Real lam){
  register int i; /*color*/
  register int j;
  register site* s;
  register complex ctmp0,ctmp1;
  const Real lam05=lam*0.5;
  const Real lamsr05=lam*SR05;

#ifdef PREFETCH
  const int boundfetch=sites_on_node-FETCH_UP;
#endif
  half_wilson_vector sd;
  
  su3_matrix * link=&t_blocked_link2[sites_on_node*dir];

  if(sign == PLUS){
      switch(dir){
	  case 3: /* x */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		  if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);

		  for(i=0;i<3;i++){
		      (dest[j].d[3].c[i]).real += lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[3].c[i]).imag -= lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[2].c[i]).real += lam*(sd.h[1].c[i]).imag; 
		      (dest[j].d[2].c[i]).imag -= lam*(sd.h[1].c[i]).real; 
		  }
	      }
	      break;

	  case 9: /* y */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		  if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){
		      (dest[j].d[3].c[i]).real -= lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[3].c[i]).imag -= lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[2].c[i]).real += lam*(sd.h[1].c[i]).real; 
		      (dest[j].d[2].c[i]).imag += lam*(sd.h[1].c[i]).imag; 
		  }



	      }
	      break;

	  case 13: /* z */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);

		  for(i=0;i<3;i++){
		      (dest[j].d[2].c[i]).real += lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[2].c[i]).imag -= lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[3].c[i]).real -= lam*(sd.h[1].c[i]).imag; 
		      (dest[j].d[3].c[i]).imag += lam*(sd.h[1].c[i]).real; 
		  }

	      }
	      break;

	  case 15: /* t */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  scalar_mult_add_su3_vector(&dest[j].d[2],&sd.h[0],lam,&dest[j].d[2]);
		  scalar_mult_add_su3_vector(&dest[j].d[3],&sd.h[1],lam,&dest[j].d[3]);
	      }
	      break;

	      /* diagonal offsets */
	  case 0: /*  1 1 0 0 */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){
		      dest[j].d[2].c[i].real += lamsr05*( sd.h[1].c[i].real  +sd.h[1].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( -sd.h[1].c[i].real  +sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].real += lamsr05*( -sd.h[0].c[i].real  +sd.h[0].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( -sd.h[0].c[i].real  -sd.h[0].c[i].imag);
		  }
	      }
	      break;

	  case  1:/*  1 0 1 0 */

	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){

		      ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		      ctmp1.imag = lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real += SR05*(ctmp0.imag + ctmp1.imag);
		      dest[j].d[2].c[i].imag += -SR05*(ctmp0.real + ctmp1.real);

		      dest[j].d[3].c[i].real += -SR05*(ctmp1.imag - ctmp0.imag);
		      dest[j].d[3].c[i].imag += SR05*(ctmp1.real - ctmp0.real);
		  }
	      }
	      break;

	  case 2: /* 1 0 0 1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real = lam05*(sd.h[0].c[i].real + sd.h[1].c[i].real);
		      ctmp0.imag = lam05*(sd.h[0].c[i].imag + sd.h[1].c[i].imag);
		      ctmp1.real = lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);
		      
		      dest[j].d[2].c[i].real +=  SR05*(ctmp0.real + ctmp1.imag);
		      dest[j].d[2].c[i].imag +=  SR05*(ctmp0.imag - ctmp1.real);
		      dest[j].d[3].c[i].real +=  SR05*(ctmp1.real + ctmp0.imag);
		      dest[j].d[3].c[i].imag +=  SR05*(ctmp1.imag - ctmp0.real);

		  }
	      }
	      break;

	  case 4: /* 1 0 0 -1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real =lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real += SR05*(-ctmp0.real + ctmp1.imag);
		      dest[j].d[2].c[i].imag += - SR05*(ctmp0.imag + ctmp1.real);

		      dest[j].d[3].c[i].real += SR05*(-ctmp1.real + ctmp0.imag);
		      dest[j].d[3].c[i].imag += SR05*(-ctmp1.imag - ctmp0.real);

		  }
	      }
	      break;

	  case 5: /* 1 0 -1 0 */


	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      
		      ctmp0.real=  lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag=  lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real=  lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		      ctmp1.imag=  lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real += SR05*(ctmp1.imag - ctmp0.imag);
		      dest[j].d[2].c[i].imag += -SR05*(ctmp1.real - ctmp0.real);

		      dest[j].d[3].c[i].real += SR05*(ctmp1.imag + ctmp0.imag);
		      dest[j].d[3].c[i].imag += -SR05*(ctmp1.real + ctmp0.real);
		  }
	      }
	      break;

	  case 6: /* 1 -1 0 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( -sd.h[1].c[i].real  +sd.h[1].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( -sd.h[1].c[i].real  -sd.h[1].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( sd.h[0].c[i].real  +sd.h[0].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( -sd.h[0].c[i].real  +sd.h[0].c[i].imag);

		  }
	      }
	      break;

	  case 7:/*  0 1 1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real = lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag = lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real = lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);


		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real+= SR05*(ctmp1.real + ctmp0.imag);
		      dest[j].d[2].c[i].imag += SR05*(ctmp1.imag - ctmp0.real);

		      dest[j].d[3].c[i].real+= SR05*(-ctmp0.real - ctmp1.imag);
		      dest[j].d[3].c[i].imag += SR05*(-ctmp0.imag + ctmp1.real);

		  }
	      }

	      break;


	  case 8: /* 0 1 0 1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){

		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		  ctmp1.imag = lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);

		  dest[j].d[2].c[i].real += SR05*(ctmp1.real + ctmp0.real);
		  dest[j].d[2].c[i].imag += SR05*(ctmp1.imag + ctmp0.imag);

		  dest[j].d[3].c[i].real += SR05*(ctmp1.real - ctmp0.real);
		  dest[j].d[3].c[i].imag += SR05*(ctmp1.imag - ctmp0.imag);

	      }
	      }
	      break;


	  case 10: /* 0 1 0 -1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){
		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		  ctmp1.imag =lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);

		  dest[j].d[2].c[i].real+= SR05*(ctmp1.real - ctmp0.real);
		  dest[j].d[2].c[i].imag += SR05*(ctmp1.imag - ctmp0.imag);

		  dest[j].d[3].c[i].real+= -SR05*(ctmp1.real + ctmp0.real);
		  dest[j].d[3].c[i].imag += -SR05*(ctmp1.imag + ctmp0.imag);
	      }
	      }
	      break;


	  case 11: /* 0 1 -1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){

		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real =lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		  ctmp1.imag =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);
		  
		  dest[j].d[2].c[i].real += SR05*(ctmp1.real - ctmp0.imag);
		  dest[j].d[2].c[i].imag += SR05*(ctmp1.imag + ctmp0.real);

		  dest[j].d[3].c[i].real +=  SR05*(-ctmp0.real + ctmp1.imag);
		  dest[j].d[3].c[i].imag += SR05*(-ctmp0.imag - ctmp1.real);
	      }
	      }
	      break;


	  case 12: /* 0 0 1 1 */

	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( sd.h[0].c[i].real  +sd.h[0].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( -sd.h[0].c[i].real  +sd.h[0].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( sd.h[1].c[i].real  -sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( sd.h[1].c[i].real  +sd.h[1].c[i].imag);

		  }
	      }
	      break;

	  case 14: /* 0 0 1 -1 */


	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( -sd.h[0].c[i].real  +sd.h[0].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( -sd.h[0].c[i].real  -sd.h[0].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( -sd.h[1].c[i].real  -sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( sd.h[1].c[i].real  -sd.h[1].c[i].imag);

		  }
	      }
	      break;


	  default:
	      node0_printf("BAD CALL TO WP_GROW()\n");
      }

  }
  if(sign == MINUS){
      switch(dir){

	  case 3: /* XDOWN: */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);

		  for(i=0;i<3;i++){
		      (dest[j].d[3].c[i]).real -= lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[3].c[i]).imag += lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[2].c[i]).real -= lam*(sd.h[1].c[i]).imag; 
		      (dest[j].d[2].c[i]).imag += lam*(sd.h[1].c[i]).real; 
		  }
	      }

	      break;

	  case 9: /* YDOWN: */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){
		      (dest[j].d[3].c[i]).real += lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[3].c[i]).imag += lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[2].c[i]).real -= lam*(sd.h[1].c[i]).real; 
		      (dest[j].d[2].c[i]).imag -= lam*(sd.h[1].c[i]).imag; 
		  }
	      }

	      break;

	  case 13: /* ZDOWN: */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);

		  for(i=0;i<3;i++){
		      (dest[j].d[2].c[i]).real -= lam*(sd.h[0].c[i]).imag; 
		      (dest[j].d[2].c[i]).imag += lam*(sd.h[0].c[i]).real; 
		      (dest[j].d[3].c[i]).real += lam*(sd.h[1].c[i]).imag; 
		      (dest[j].d[3].c[i]).imag -= lam*(sd.h[1].c[i]).real; 
		  }

	      }

	      break;

	  case 15: /* TDOWN: */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  scalar_mult_add_su3_vector(&dest[j].d[2],&sd.h[0],-lam,&dest[j].d[2]);
		  scalar_mult_add_su3_vector(&dest[j].d[3],&sd.h[1],-lam,&dest[j].d[3]);
	      }

	      break;

	      /* diagonal offsets */
	  case 0: /*  1 1 0 0 */
	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){
		      dest[j].d[2].c[i].real += lamsr05*( -sd.h[1].c[i].real  -sd.h[1].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( sd.h[1].c[i].real  -sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].real += lamsr05*( sd.h[0].c[i].real  -sd.h[0].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( sd.h[0].c[i].real  +sd.h[0].c[i].imag);
		  }
	      }
	      break;

	  case  1:/*  1 0 1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){

		      ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		      ctmp1.imag = lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real += -SR05*(ctmp0.imag + ctmp1.imag);
		      dest[j].d[2].c[i].imag += SR05*(ctmp0.real + ctmp1.real);

		      dest[j].d[3].c[i].real += SR05*(ctmp1.imag - ctmp0.imag);
		      dest[j].d[3].c[i].imag += -SR05*(ctmp1.real - ctmp0.real);
		  }
	      }
		  break;

	  case 2: /* 1 0 0 1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real = lam05*(sd.h[0].c[i].real + sd.h[1].c[i].real);
		      ctmp0.imag = lam05*(sd.h[0].c[i].imag + sd.h[1].c[i].imag);
		      ctmp1.real = lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);
		      
		      dest[j].d[2].c[i].real -=  SR05*(ctmp0.real + ctmp1.imag);
		      dest[j].d[2].c[i].imag -=  SR05*(ctmp0.imag - ctmp1.real);
		      dest[j].d[3].c[i].real -=  SR05*(ctmp1.real + ctmp0.imag);
		      dest[j].d[3].c[i].imag -=  SR05*(ctmp1.imag - ctmp0.real);

		  }
	      }
		  break;

	  case 4: /* 1 0 0 -1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real =lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real -= SR05*(-ctmp0.real + ctmp1.imag);
		      dest[j].d[2].c[i].imag += SR05*(ctmp0.imag + ctmp1.real);

		      dest[j].d[3].c[i].real -= SR05*(-ctmp1.real + ctmp0.imag);
		      dest[j].d[3].c[i].imag += SR05*(+ctmp1.imag + ctmp0.real);

		  }
	      }
	      break;
	  case 5: /* 1 0 -1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      
		      ctmp0.real=  lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag=  lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real=  lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		      ctmp1.imag=  lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real -= SR05*(ctmp1.imag - ctmp0.imag);
		      dest[j].d[2].c[i].imag += SR05*(ctmp1.real - ctmp0.real);

		      dest[j].d[3].c[i].real -= SR05*(ctmp1.imag + ctmp0.imag);
		      dest[j].d[3].c[i].imag += SR05*(ctmp1.real + ctmp0.real);
		  }
	      }
	      break;

	  case 6: /* 1 -1 0 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( sd.h[1].c[i].real - sd.h[1].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( sd.h[1].c[i].real + sd.h[1].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( -sd.h[0].c[i].real  -sd.h[0].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( sd.h[0].c[i].real  -sd.h[0].c[i].imag);

		  }
	      }
		  break;

	  case 7:/*  0 1 1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  for(i=0;i<3;i++){
		      ctmp0.real = lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		      ctmp0.imag = lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		      ctmp1.real = lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		      ctmp1.imag = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		      CSUM(dest[j].d[0].c[i],ctmp0);
		      CSUM(dest[j].d[1].c[i],ctmp1);

		      dest[j].d[2].c[i].real -= SR05*(ctmp1.real + ctmp0.imag);
		      dest[j].d[2].c[i].imag -= SR05*( ctmp1.imag - ctmp0.real);

		      dest[j].d[3].c[i].real += SR05*( ctmp0.real + ctmp1.imag);
		      dest[j].d[3].c[i].imag -= SR05*(-ctmp0.imag + ctmp1.real);

		  }
	      }
		  break;

	  case 8: /* 0 1 0 1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){

		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real = lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		  ctmp1.imag = lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);

		  dest[j].d[2].c[i].real -= SR05*(ctmp1.real + ctmp0.real);
		  dest[j].d[2].c[i].imag -= SR05*(ctmp1.imag + ctmp0.imag);

		  dest[j].d[3].c[i].real -= SR05*(ctmp1.real - ctmp0.real);
		  dest[j].d[3].c[i].imag -= SR05*(ctmp1.imag - ctmp0.imag);

	      }
	      }
	      break;

	  case 10: /* 0 1 0 -1 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){
		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);
		  ctmp1.imag =lam05*(-sd.h[0].c[i].real + sd.h[1].c[i].real);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);

		  dest[j].d[2].c[i].real -= SR05*(ctmp1.real - ctmp0.real);
		  dest[j].d[2].c[i].imag -= SR05*(ctmp1.imag - ctmp0.imag);

		  dest[j].d[3].c[i].real += SR05*(ctmp1.real + ctmp0.real);
		  dest[j].d[3].c[i].imag += SR05*(ctmp1.imag + ctmp0.imag);
	      }
	      }
	      break;

	  case 11: /* 0 1 -1 0 */

	      FORALLSITES(j,s)
	      {
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
	      for(i=0;i<3;i++){

		  ctmp0.real =lam05*(sd.h[0].c[i].real+sd.h[1].c[i].real);
		  ctmp0.imag =lam05*(sd.h[0].c[i].imag+sd.h[1].c[i].imag);

		  ctmp1.real =lam05*(sd.h[0].c[i].real - sd.h[1].c[i].real);
		  ctmp1.imag =lam05*(sd.h[0].c[i].imag - sd.h[1].c[i].imag);

		  CSUM(dest[j].d[0].c[i],ctmp0);
		  CSUM(dest[j].d[1].c[i],ctmp1);
		  
		  dest[j].d[2].c[i].real-= SR05*(ctmp1.real - ctmp0.imag);
		  dest[j].d[2].c[i].imag-= SR05*(ctmp1.imag + ctmp0.real);

		  dest[j].d[3].c[i].real-= SR05*(-ctmp0.real + ctmp1.imag);
		  dest[j].d[3].c[i].imag+= SR05*( ctmp0.imag + ctmp1.real);
	      }
	      }
	      break;

	  case 12: /* 0 0 1 1 */

	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( -sd.h[0].c[i].real  -sd.h[0].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*(  sd.h[0].c[i].real  -sd.h[0].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( -sd.h[1].c[i].real  +sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( -sd.h[1].c[i].real  -sd.h[1].c[i].imag);

		  }
	      }
	      break;
	  case 14: /* 0 0 1 -1 */


	      FORALLSITES(j,s){
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_WW(&src[j+FETCH_UP],&dest[j+FETCH_UP]);}
#endif
	          mult_su3_mat_hwvec( &link[j], (half_wilson_vector * )(gen_pt[0][j]), &sd ); 
		  scalar_mult_add_su3_vector(&dest[j].d[0],&sd.h[0],lam,&dest[j].d[0]);
		  scalar_mult_add_su3_vector(&dest[j].d[1],&sd.h[1],lam,&dest[j].d[1]);
		  for(i=0;i<3;i++){

		      dest[j].d[2].c[i].real += lamsr05*( sd.h[0].c[i].real  -sd.h[0].c[i].imag);
		      dest[j].d[2].c[i].imag += lamsr05*( sd.h[0].c[i].real  +sd.h[0].c[i].imag);

		      dest[j].d[3].c[i].real += lamsr05*( sd.h[1].c[i].real  +sd.h[1].c[i].imag);
		      dest[j].d[3].c[i].imag += lamsr05*( -sd.h[1].c[i].real  +sd.h[1].c[i].imag);

		  }
	      }
	      break;
	  default:
		  printf("BAD CALL TO WP_GROW()\n");
      }
  }
}

