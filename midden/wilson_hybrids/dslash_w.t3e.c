/************* dslash_w.t3e.c *******************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

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
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i,end;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];
register int first,count;

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
      if( i < loopend-FETCH_UP ){
	prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
      }
      wp_shrink_4dir( F_PT(s,src), &(s->htmp[XUP]), &(s->htmp[YUP]),
		      &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    } END_LOOP
    for( dir=XUP; dir <= TUP; dir++) {
	tag[dir]=start_gather( F_OFFSET(htmp[dir]), sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
      if( i < loopend-FETCH_UP ){
	prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
	prefetch_M( &((s+FETCH_UP)->link[XUP]) );
	prefetch_M( &((s+FETCH_UP)->link[YUP]) );
	prefetch_M( &((s+FETCH_UP)->link[ZUP]) );
	prefetch_M( &((s+FETCH_UP)->link[TUP]) );
      }
      wp_shrink_4dir( F_PT(s,src), &hwvx, &hwvy, &hwvz, &hwvt, -isign);
      mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
      mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
      mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
      mult_adj_su3_mat_hwvec( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    } END_LOOP

    for( dir=XUP; dir <= TUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather(F_OFFSET(htmp[OPP_DIR(dir)]), 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_4MWWWW( 
			&((s+FETCH_UP)->link[XUP]), 
			(half_wilson_vector *)(gen_pt[XUP][i+FETCH_UP]),
			(half_wilson_vector *)(gen_pt[YUP][i+FETCH_UP]),
			(half_wilson_vector *)(gen_pt[ZUP][i+FETCH_UP]),
			(half_wilson_vector *)(gen_pt[TUP][i+FETCH_UP]) );
      }
      mult_su3_mat_hwvec( &(s->link[XUP]), 
			  (half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
      mult_su3_mat_hwvec( &(s->link[YUP]), 
			  (half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
      mult_su3_mat_hwvec( &(s->link[ZUP]), 
			  (half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
      mult_su3_mat_hwvec( &(s->link[TUP]), 
			  (half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
      grow_add_four_wvecs( F_PT(s,dest),
			   &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    } END_LOOP

    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_W( (half_wilson_vector *)(gen_pt[XDOWN][i+FETCH_UP]) );
	prefetch_W( (half_wilson_vector *)(gen_pt[YDOWN][i+FETCH_UP]) );
	prefetch_W( (half_wilson_vector *)(gen_pt[ZDOWN][i+FETCH_UP]) );
	prefetch_W( (half_wilson_vector *)(gen_pt[TDOWN][i+FETCH_UP]) );
	prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
      } 
      grow_add_four_wvecs( F_PT(s,dest),
			   (half_wilson_vector *)(gen_pt[XDOWN][i]),
			   (half_wilson_vector *)(gen_pt[YDOWN][i]),
			   (half_wilson_vector *)(gen_pt[ZDOWN][i]),
			   (half_wilson_vector *)(gen_pt[TDOWN][i]),
			   -isign, 1 );	/* "1" SUMs in current dest */
    } END_LOOP

    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash() ) */


/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=1,use
  start_gather, otherwise use restart_gather.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_special(field_offset src,field_offset dest,int isign,
		    int parity,msg_tag **tag,int is_started)
{
wilson_vector localvec;	/* temporary storage */
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;
register int first,count;

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir( F_PT(s,src), &(s->htmp[XUP]), &(s->htmp[YUP]),
            &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    } END_LOOP
    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[dir]=start_gather( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir] );
	else restart_gather( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir], tag[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir( F_PT(s,src), &hwvx, &hwvy, &hwvz, &hwvt, -isign);
	mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    } END_LOOP

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[OPP_DIR(dir)]=start_gather(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	mult_su3_mat_hwvec( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	mult_su3_mat_hwvec( &(s->link[TUP]), 
		(half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
	grow_add_four_wvecs( F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    } END_LOOP

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	grow_add_four_wvecs( F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    } END_LOOP

} /* end (of dslash_special() ) */
