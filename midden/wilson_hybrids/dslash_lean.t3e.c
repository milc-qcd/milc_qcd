/************* dslash_lean.t3e.c *******************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Version to use less memory - do things one direction at a time */

/*
	dslash(F_OFFSET(psi),F_OFFSET(mp),isign,l_parity);
Compute SUM_dirs ( 
    ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
)

*/

#include "wi_hyb_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

void dslash(field_offset src,field_offset dest,int isign,int parity)
{
wilson_vector localvec;	/* temporary storage */
half_wilson_vector hwv;

register int i,end;
register site *s;
register int dir,otherparity;
msg_tag *tag[2];
register int first,count;

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    for(dir=XUP;dir<=TUP;dir++){
        /* Take Wilson projection for src displaced in up direction, gather
           it to "our site" */
        FORSOMEPARITY(i,s,otherparity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
	  }
	  wp_shrink( (wilson_vector *)F_PT(s,src), 
		     &(s->htmp[0]), dir, isign);
        } END_LOOP
	tag[0]=start_gather( F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
	    dir, parity, gen_pt[0] );

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
        FORSOMEPARITY(i,s,otherparity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
	    prefetch_M( &((s+FETCH_UP)->link[dir]) );
	  }
	  wp_shrink( (wilson_vector *)F_PT(s,src), 
		     &hwv, dir, -isign);
	  mult_adj_su3_mat_hwvec( &(s->link[dir]), &hwv, &(s->htmp[1]));
        } END_LOOP

	tag[1]=start_gather(F_OFFSET(htmp[1]), 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[1] );


        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
	wait_gather(tag[0]);

        if(dir==XUP)FORSOMEPARITY(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	    prefetch_H( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i+FETCH_UP]) );
	  }
	  mult_su3_mat_hwvec( &(s->link[dir]), 
			      (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	  wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
        } END_LOOP
        else FORSOMEPARITY(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	    prefetch_H( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i+FETCH_UP]) );
	  }
	  mult_su3_mat_hwvec( &(s->link[dir]), 
			      (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	  wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
        } END_LOOP
	cleanup_gather(tag[0]);

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
	wait_gather(tag[1]);

        FORSOMEPARITY(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	    prefetch_H( (half_wilson_vector *)(gen_pt[1][i+FETCH_UP]) );
	  }
	  wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		       (wilson_vector *)F_PT(s,dest), dir, -isign);
        } END_LOOP
	cleanup_gather(tag[1]);

    } /* end loop over directions */
} /* end (of dslash_w.c) */
