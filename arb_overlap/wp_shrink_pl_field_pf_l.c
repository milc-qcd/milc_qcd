/************* wp_shrink_pl_field_pf_l.c  **************************/
/* 
   Compute the "Wilson projection" of a Wilson fermion vector.
   (1 +- n dot gamma_j) is a projection operator, and we want to isolate
   the components of the vector that it keeps.  In other words, keep
   the components of the vector along the eigenvectors of 1+-gamma_j
   with eigenvalue 2, and throw away those with eigenvalue 0.

usage:  wp_shrink_pl( wilson_vector *src, half_wilson_vector *dest,
int offset, int sign )



Planar offset table

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
#include <arb_ov_includes.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif


#define SR05 0.707106781187
#define SR08 0.353553390593


#ifdef PREFETCH
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#endif



void wp_shrink_pl_field_l( wilson_vector  * __restrict__ src, half_wilson_vector  * __restrict__ dest,
	const int dir, const int sign ){
    register int i; /*color*/
    register int j;
    register site* s;
    complex tvec2,tvec4;
    Real sr05;
    su3_vector vtvec1,vtvec2,vtvec3,vtvec5;
    half_wilson_vector dd;
#ifdef PREFETCH
    const int  boundfetch=sites_on_node-FETCH_UP;
#endif
    su3_matrix * link=&t_blocked_link2[sites_on_node*dir];
    if(sign == PLUS){
	switch(dir){
	    case 3: /* x */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real - src[j].d[3].c[i].imag;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag + src[j].d[3].c[i].real;
			dd.h[1].c[i].real = src[j].d[1].c[i].real - src[j].d[2].c[i].imag;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag + src[j].d[2].c[i].real;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;
	    case  9: /* y */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real - src[j].d[3].c[i].real;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag - src[j].d[3].c[i].imag;
			dd.h[1].c[i].real = src[j].d[1].c[i].real + src[j].d[2].c[i].real;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag + src[j].d[2].c[i].imag;
		    }   
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;
	    case 13: /* z */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real - src[j].d[2].c[i].imag;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag + src[j].d[2].c[i].real;
			dd.h[1].c[i].real = src[j].d[1].c[i].real + src[j].d[3].c[i].imag;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag - src[j].d[3].c[i].real;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;
	    case 15: /* t */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real + src[j].d[2].c[i].real;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag + src[j].d[2].c[i].imag;
			dd.h[1].c[i].real = src[j].d[1].c[i].real + src[j].d[3].c[i].real;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag + src[j].d[3].c[i].imag;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

		/* diagonal offsets */
	    case 0: /*  1 1 0 0 */

		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real= ( -src[j].d[3].c[i].real  -src[j].d[3].c[i].imag);
			vtvec1.c[i].imag= ( src[j].d[3].c[i].real  -src[j].d[3].c[i].imag);
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]) );

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[2].c[i].real  + src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case  1:/*  1 0 1 0 */

		sr05=SR05;
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  -tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));


		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real=  tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  +tvec4.imag ;

		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}

		break;

	    case 2: /* 1 0 0 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  tvec2.real  -tvec2.imag  ;
			vtvec3.c[i].imag=  tvec2.real  +tvec2.imag  ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;

			vtvec3.c[i].real=        tvec2.real  +tvec2.imag  ;
			vtvec3.c[i].imag=        -tvec2.real  +tvec2.imag  ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 4: /* 1 0 0 -1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}
		break;

	    case 5: /* 1 0 -1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  -tvec4.real  +tvec4.imag;
			vtvec5.c[i].imag=  -tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real=  tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 6: /* 1 -1 0 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  -src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[2].c[i].real  - src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 7:/*  0 1 1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 8: /* 0 1 0 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;


			vtvec5.c[i].real=  tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 10: /* 0 1 0 -1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real= -tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &(vtvec2),&vtvec5,sr05,&(dd.h[0])  );

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  -tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 11: /* 0 1 -1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;
			vtvec3.c[i].real=  tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;




	    case 12: /* 0 0 1 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
			vtvec1.c[i].imag= src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 14: /* 0 0 1 -1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  -src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  -src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;



	    default:
		node0_printf("BAD CALL TO WP_SHRINK()\n");
	}
    }

    if(sign == MINUS){
	switch(dir){
	    case 3: /* XDOWN: */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real + src[j].d[3].c[i].imag;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag - src[j].d[3].c[i].real;
			dd.h[1].c[i].real = src[j].d[1].c[i].real + src[j].d[2].c[i].imag;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag - src[j].d[2].c[i].real;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 9: /*YDOWN: */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real + src[j].d[3].c[i].real;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag + src[j].d[3].c[i].imag;
			dd.h[1].c[i].real = src[j].d[1].c[i].real - src[j].d[2].c[i].real;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag - src[j].d[2].c[i].imag;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 13: /*ZDOWN: */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real + src[j].d[2].c[i].imag;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag - src[j].d[2].c[i].real;
			dd.h[1].c[i].real = src[j].d[1].c[i].real - src[j].d[3].c[i].imag;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag + src[j].d[3].c[i].real;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 15: /*TDOWN: */
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			dd.h[0].c[i].real = src[j].d[0].c[i].real - src[j].d[2].c[i].real;
			dd.h[0].c[i].imag = src[j].d[0].c[i].imag - src[j].d[2].c[i].imag;
			dd.h[1].c[i].real = src[j].d[1].c[i].real - src[j].d[3].c[i].real;
			dd.h[1].c[i].imag = src[j].d[1].c[i].imag - src[j].d[3].c[i].imag;
		    }
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;


		/* diagonal offsets */

	    case 0: /*  1 1 0 0 */

		sr05=SR05;
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real= ( src[j].d[3].c[i].real  +src[j].d[3].c[i].imag);
			vtvec1.c[i].imag= ( -src[j].d[3].c[i].real  +src[j].d[3].c[i].imag);
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real= ( -src[j].d[2].c[i].real  +src[j].d[2].c[i].imag);
			vtvec1.c[i].imag= ( -src[j].d[2].c[i].real  - src[j].d[2].c[i].imag);
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;


	    case  1:/*  1 0 1 0 */

		sr05=SR05;
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real=  -tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 2: /* 1 0 0 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  -tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 4: /* 1 0 0 -1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;
			vtvec3.c[i].real= tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;


	    case 5: /* 1 0 -1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;

			vtvec5.c[i].real=  tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;
			vtvec5.c[i].real=  -tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 6: /* 1 -1 0 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){

			vtvec1.c[i].real=  -src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}

		break;

	    case 7:/*  0 1 1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;
			vtvec3.c[i].real= -tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  -tvec2.imag;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  tvec2.real  +tvec2.imag ;
			vtvec3.c[i].imag=  -tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}
		break;


	    case 8: /* 0 1 0 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real=  -tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;
			vtvec5.c[i].real=  -tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  -tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}
		break;


	    case 10: /* 0 1 0 -1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real - src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].real;

			vtvec5.c[i].real= tvec4.real  +tvec4.imag ;
			vtvec5.c[i].imag=  -tvec4.real  +tvec4.imag;
		    }
		    scalar_mult_add_su3_vector( &(vtvec2),&vtvec5,sr05,&(dd.h[0])  );

		    for(i=0;i<3;i++){
			vtvec2.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].imag;
			vtvec2.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].real;

			tvec4.real = src[j].d[2].c[i].real + src[j].d[3].c[i].imag;
			tvec4.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].real;
			vtvec5.c[i].real=  tvec4.real  -tvec4.imag ;
			vtvec5.c[i].imag=  tvec4.real  +tvec4.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec2,&vtvec5,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}
		break;


	    case 11: /* 0 1 -1 0 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real + src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag + src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real - src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag - src[j].d[3].c[i].imag;
			vtvec3.c[i].real=  -tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  -tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real = src[j].d[0].c[i].real - src[j].d[1].c[i].real;
			vtvec1.c[i].imag = src[j].d[0].c[i].imag - src[j].d[1].c[i].imag;

			tvec2.real = src[j].d[2].c[i].real + src[j].d[3].c[i].real;
			tvec2.imag = src[j].d[2].c[i].imag + src[j].d[3].c[i].imag;

			vtvec3.c[i].real=  tvec2.real  -tvec2.imag ;
			vtvec3.c[i].imag=  tvec2.real  +tvec2.imag ;
		    }
		    scalar_mult_add_su3_vector( &vtvec1,&vtvec3,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);

		}
		break;


	    case 12: /* 0 0 1 1 */
		sr05=SR05;

		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  -src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[2].c[i].real  -src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  -src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;

	    case 14: /* 0 0 1 -1 */

		sr05=SR05;
		FORALLSITES(j,s)
		{
#ifdef PREFETCH
		    if(j<boundfetch) {prefetch_M(&link[j+FETCH_UP]);prefetch_W(&src[j+FETCH_UP]);}
#endif
		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
			vtvec1.c[i].imag=  -src[j].d[2].c[i].real  +src[j].d[2].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[0]),&vtvec1,sr05,&(dd.h[0]));

		    for(i=0;i<3;i++){
			vtvec1.c[i].real=  src[j].d[3].c[i].real  -src[j].d[3].c[i].imag;
			vtvec1.c[i].imag=  src[j].d[3].c[i].real  +src[j].d[3].c[i].imag;
		    }
		    scalar_mult_add_su3_vector( &(src[j].d[1]),&vtvec1,sr05,&(dd.h[1]));
		    mult_adj_su3_mat_hwvec(&link[j],&dd,&dest[j]);
		}
		break;


	    default:
		printf("BAD CALL TO WP_SHRINK()\n");
	}
    }

}

