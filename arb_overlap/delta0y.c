/*************** delta0y.c ***********************/
/* MIMD version 7 */

/*!\file delta0y.c
 *\brief multiply the (kernel) Dirac operator onto a vector
 */

/*! \fn void delta0_field (wilson_vector *  src,wilson_vector *  dest,  int isign)
 * \brief multiply the (kernel) Dirac operator onto a vector
 * \param src vector on which the Dirac operator is to be multiplied
 * \param dest result d  src
 * \param isign isign=1 for \f$ d_0 \f$, -1 for \f$ d_0^\dagger \f$
 *  This version uses gathers, not general gathers--
 *  the gather offsets are initialized in setup_offsets.c
 */

/* MIMD version 6 */

/*
   dest= delta0*src
*/



#include "arb_ov_includes.h"

/* prefetching for field-wise...*/

#ifdef PREFETCH
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#endif


#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif


void wp_shrink_pl_field( wilson_vector *src, half_wilson_vector *dest,
	int dir, int sign);


void wp_shrink_pl_field_l( wilson_vector *src, half_wilson_vector *dest,
	int dir, int sign);


void wp_grow_pl_field( half_wilson_vector *dest, wilson_vector *src,
	int dir, int sign  , Real lam);
void wp_grow_pl_field_l( half_wilson_vector *dest, wilson_vector *src,
	int dir, int sign  , Real lam);

void delta0_field(wilson_vector *  src,wilson_vector *  dest,  int isign)
{
    register int i;
    register site *s;
    msg_tag *tag[2];
    int ipath,k;
    half_wilson_vector *  htmp1, *htmp2, *htmp3;
    wilson_vector *  tmp1;

    htmp1=delta0_htmp1;
    htmp2=delta0_htmp2;
    htmp3=delta0_htmp3;
    tmp1=delta0_tmp1;

    ndelta0++;

    if(clover_term == 0.0){
	/* begin with the origin */
	FORALLSITES(i,s) {
#ifdef PREFETCH
	    if(i< sites_on_node - FETCH_UP){
		prefetch_W(&(src[i+FETCH_UP]));
		prefetch_W(&(dest[i+FETCH_UP]));
	    }
#endif
	    scalar_mult_wvec(&(src[i]), lambda[0],&(dest[i]));
	}
    } else {
	
    /* multiply the clover term, if needed */
	mult_ldu1_field(src, tmp1, t_clov,t_clov_diag, EVENANDODD);
	
	FORALLSITES(i,s) {
#ifdef PREFETCH
	    if(i< sites_on_node - FETCH_UP){
		prefetch_W(&(src[i+FETCH_UP]));
		prefetch_W(&(tmp1[i+FETCH_UP]));
		prefetch_W(&(dest[i+FETCH_UP]));
	    }
#endif
	    scalar_mult_add_wvec(&(tmp1[i]),&(src[i]),lambda[0],&(dest[i]));
	}
    }
    /* now loop over the offsets and collect them */
    for(ipath=0;ipath<off_max;ipath++){
	k=label[ipath];
	/*node0_printf("ipath=%d k=%d\n",ipath,k); */

	/*
	   printf("ipath=%d k=%d  %e   %e \n",ipath,k,rho[k],lambda[k]);
	   */
	/* shrink source */
	wp_shrink_pl_field(src, htmp1, ipath,isign);

	/* Take  src displaced in up direction, gather it to "our site" */

	tag[0] = start_gather_field( htmp1,
		sizeof(half_wilson_vector), goffset[ipath],
		EVENANDODD, gen_pt[0] );

	/* Take  src displaced in down direction,
	   multiply it by adjoint link matrix, gather it "up" */

	wp_shrink_pl_field_l(src, htmp2, ipath,-isign);
	tag[1] = start_gather_field( htmp2,
		sizeof(half_wilson_vector), goffset[ipath]+1,
		EVENANDODD, gen_pt[1] );



	/* Take  src displaced in up direction, gathered,
	   multiply it by link matrix,  and add to dest */
	wait_gather(tag[0]);
	/*node0_printf("mult 1\n"); */
	wp_grow_pl_field_l( htmp3, dest, ipath, isign,lambda[k]);
	cleanup_gather(tag[0]);
	
	/* multiply by the term behind*/
	wait_gather(tag[1]);

	FORALLSITES(i,s) htmp3[i]=*(half_wilson_vector *)gen_pt[1][i];
	wp_grow_pl_field( htmp3, dest, ipath, -isign,lambda[k]);
	cleanup_gather(tag[1]);
    } /* offsets */

}
/* delta0_field.c */


/* should be in libraries */

/********************  hw_to_w.c  (in su3.a) ********************
*
void hw_to_w(half_wilson_vector *src1, int ichiral,
        wilson_vector *dest) {

*Build a chiral Wilson vector from a half_wilson_vector
* dest  <-  src1
*/

void hw_to_w(half_wilson_vector *src1, int ichiral,
        wilson_vector *dest) {

  register int i,j,oppchiral;
  for(i=0;i<2;i++)for(j=0;j<3;j++)
        dest->d[i+ichiral].c[j] = src1->h[i].c[j];

  oppchiral= (2+ichiral) %4 ;
  for(i=0;i<2;i++)for(j=0;j<3;j++)
    dest->d[i+oppchiral].c[j] = cmplx(0.0,0.0);



}
/********************  w_to_hw.c  (in su3.a) ********************
*
void w_to_hw(wilson_vector *src1, int ichiral,
        half_wilson_vector *dest) {

*  Take one chirality out of a Wilson vector and make it a half_wilson_vector
* dest  <-  src1
*/


void w_to_hw(wilson_vector *src1, int ichiral, half_wilson_vector *dest)
{

  register int i;


  for(i=0;i<2;i++){
    dest->h[i]  = src1->d[i+ichiral];
  }

}


Real magsq_hwvec( half_wilson_vector *vec ){
  register int i;
  register Real sum;
  sum=0.0;
  for(i=0;i<2;i++)sum += magsq_su3vec( &(vec->h[i]) );
  return(sum);


}
