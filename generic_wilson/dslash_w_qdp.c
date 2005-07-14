/*************  dslash_w_site.c *******************************/
/* MIMD version 7 */

/*
	dslash_w_site(F_OFFSET(psi),F_OFFSET(mp),isign,l_parity);
Compute SUM_dirs ( 
    ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
)

*/

#include "generic_wilson_includes.h"
#include <qdp.h>
#include "../include/generic_qdp.h"

//#define MILC_GAMMA

static void
my_wp_shrink( wilson_vector *src, half_wilson_vector *dest,
	      int dir, int sign )
{
#ifdef MILC_GAMMA
    wp_shrink( src, dest, dir, sign);
#else
    QLA_H_eq_spproj_D((QLA_HalfFermion *)dest, (QLA_DiracFermion *)src, dir, sign);
#endif
}

static void
my_wp_grow( half_wilson_vector *src, wilson_vector *dest,
	 int dir, int sign )
{
#ifdef MILC_GAMMA
    wp_grow( src, dest, dir, sign);
#else
    QLA_D_eq_sprecon_H((QLA_DiracFermion *)dest, (QLA_HalfFermion *)src, dir, sign);
#endif
}

static void
my_wp_grow_add( half_wilson_vector *src, wilson_vector *dest,
	     int dir, int sign )
{
#ifdef MILC_GAMMA
    wp_grow_add( src, dest, dir, sign);
#else
    QLA_D_peq_sprecon_H((QLA_DiracFermion *)dest, (QLA_HalfFermion *)src, dir, sign);
#endif
}

static void
my_wp_shrink_4dir( wilson_vector *a,  half_wilson_vector *b1,
		   half_wilson_vector *b2, half_wilson_vector *b3,
		   half_wilson_vector *b4, int sign )
{
  my_wp_shrink( a,b1,XUP,sign);
  my_wp_shrink( a,b2,YUP,sign);
  my_wp_shrink( a,b3,ZUP,sign);
  my_wp_shrink( a,b4,TUP,sign);
}

static void
my_grow_add_four_wvecs( wilson_vector *a, half_wilson_vector *b1,
		     half_wilson_vector *b2, half_wilson_vector *b3,
		     half_wilson_vector *b4, int sign, int sum )
{
  if(sum==0) my_wp_grow( b1,a,XUP,sign);
  else my_wp_grow_add( b1,a,XUP,sign);
  my_wp_grow_add( b2,a,YUP,sign);
  my_wp_grow_add( b3,a,ZUP,sign);
  my_wp_grow_add( b4,a,TUP,sign);
}

void dslash_w_site( field_offset src, field_offset dest, int isign, int parity) {
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	default:        otherparity=EVENANDODD; break;
    }

#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w_site */
    if(MAXHTMP < 8){
      printf("dslash: MAXHTMP must be 8 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
      my_wp_shrink_4dir( (wilson_vector *)F_PT(s,src), &(s->htmp[XUP]),
	  &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    }
    for( dir=XUP; dir <= TUP; dir++) {
	tag[dir]=start_gather_site( F_OFFSET(htmp[dir]), sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        my_wp_shrink_4dir( (wilson_vector *)F_PT(s,src),
	    &hwvx, &hwvy, &hwvz, &hwvt, -isign);
	mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather_site(F_OFFSET(htmp[OPP_DIR(dir)]), 
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
	mult_su3_mat_hwvec( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	mult_su3_mat_hwvec( &(s->link[TUP]), 
		(half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
	my_grow_add_four_wvecs( (wilson_vector *)F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	my_grow_add_four_wvecs( (wilson_vector *)F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w_site() ) */


/* Special dslash for use by congrad.  Uses restart_gather_site() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=0,use
  start_gather_site, otherwise use restart_gather_site.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_w_site_special(field_offset src,field_offset dest,
    int isign,int parity,msg_tag **tag,int is_started)

{
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	default:        otherparity=EVENANDODD; break;
    }

#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w_site */
    if(MAXHTMP < 8){
      printf("dslash_w_site_special: MAXHTMP must be 8 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 8){
      printf("dslash_w_site_special: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        my_wp_shrink_4dir( (wilson_vector *)F_PT(s,src), &(s->htmp[XUP]),
	    &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    }
    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[dir]=start_gather_site( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir] );
	else restart_gather_site( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir], 
			     tag[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        my_wp_shrink_4dir( (wilson_vector *)F_PT(s,src),
	    &hwvx, &hwvy, &hwvz, &hwvt, -isign);
	mult_adj_su3_mat_hwvec( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
	mult_adj_su3_mat_hwvec( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[OPP_DIR(dir)]=start_gather_site(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather_site(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)], 
	    tag[OPP_DIR(dir)] );
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
	my_grow_add_four_wvecs( (wilson_vector *)F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	my_grow_add_four_wvecs( (wilson_vector *)F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }

} /* end (of dslash_w_site_special() ) */
