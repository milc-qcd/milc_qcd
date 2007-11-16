/******* dslash.c - dslash for naive KS fermions ****/
/* MIMD version 7 */
/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/prefetch.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"

/* Temporary work space for dslash_field_special */ 
static su3_vector *temp[4] ;
static int temp_not_allocated=1 ;

static void 
cleanup_one_gather_set(msg_tag *tags[])
{
  int i;

  for(i=XUP;i<=TUP;i++){
    cleanup_gather( tags[i] );
    cleanup_gather( tags[OPP_DIR(i)] );
  }
}

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[])
{
  cleanup_one_gather_set(tags1);
  cleanup_one_gather_set(tags2);
}

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<4;i++) {
      free(temp[i]) ; 
    }
  temp_not_allocated=1 ;
}


/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_site( field_offset src, field_offset dest, int parity ) {
  msg_tag *tag[16];
  
  dslash_site_special(src, dest, parity, tag, 1 );
  cleanup_one_gather_set(tag);
}


/* Special dslash_site for use by congrad.  Uses restart_gather_site() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_site, otherwise use restart_gather_site. 
  The calling program must clean up the gathers! */
void dslash_site_special( field_offset src, field_offset dest,
			  int parity, msg_tag **tag, int start ){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  
  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Start gathers from positive directions */
  for(dir=XUP; dir<=TUP; dir++){
    if(start==1) tag[dir] = start_gather_site( src, sizeof(su3_vector),
					       dir, parity, gen_pt[dir] );
    else restart_gather_site( src, sizeof(su3_vector),
			      dir, parity, gen_pt[dir] , tag[dir] );
  }
  
  /* Multiply by adjoint matrix at other sites */
  FORSOMEPARITYDOMAIN(i,s,otherparity){
    if( i < loopend-FETCH_UP ){
      prefetch_4MV4V( &((s+FETCH_UP)->link[XUP]),
		      (su3_vector *)F_PT((s+FETCH_UP),src),
		      (s+FETCH_UP)->tempvec);
    }
    mult_adj_su3_mat_vec_4dir( s->link,
			       (su3_vector *)F_PT(s,src), s->tempvec );
  } END_LOOP
      
  /* Start gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
    if (start==1) tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
							 sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
    else restart_gather_site( F_OFFSET(tempvec[dir]), sizeof(su3_vector),
			      OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] , tag[OPP_DIR(dir)] );
  }
  
  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[dir]);
  }
  
#ifdef SCHROED_FUN
  FORSOMEPARITY(i,s,parity) if(s->t > 0){
    if(s->t == (nt-1)){
      mult_su3_mat_vec( &(s->link[XUP]),
			(su3_vector *)(gen_pt[XUP][i]), (su3_vector *)F_PT(s,dest));
      for(dir=YUP; dir<TUP; dir++){
	mult_su3_mat_vec_sum( &(s->link[dir]),
			      (su3_vector *)(gen_pt[dir][i]), (su3_vector *)F_PT(s,dest));
      }
    }
    else{
#else
      FORSOMEPARITY(i,s,parity){
	if( i < loopend-FETCH_UP ){
	  prefetch_4MVVVV(&((s+FETCH_UP)->link[XUP]), 
			  (su3_vector *)gen_pt[XUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[YUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
	}
#endif
	mult_su3_mat_vec_sum_4dir( s->link,
		   (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
		   (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
		   (su3_vector *)F_PT(s,dest));
	/*------------------------------------------------------------*/
#ifdef SCHROED_FUN
      }
#endif
  } END_LOOP
    
  /*------------------------------------------------------------*/
  /* Wait gathers from negative directions, accumulate (negative) */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[OPP_DIR(dir)]);
  }
    
  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_VVVV( 
		    (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
    }
    
    sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
		       (su3_vector *)(gen_pt[XDOWN][i]),
		       (su3_vector *)(gen_pt[YDOWN][i]),
		       (su3_vector *)(gen_pt[ZDOWN][i]),
		       (su3_vector *)(gen_pt[TDOWN][i]) ); 
      /*------------------------------------------------------------*/
  } END_LOOP
}

void dslash_field( su3_vector *src, su3_vector *dest, int parity ) {
  msg_tag *tag[16];
    
   dslash_field_special(src, dest, parity, tag, 1 );
   cleanup_one_gather_set(tag);
}


/* Special dslash_field for use by congrad.  Uses restart_gather_field() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_field, otherwise use restart_gather_field. 
  The calling program must clean up the gathers! */
void dslash_field_special( su3_vector *src, su3_vector *dest,
			  int parity, msg_tag **tag, int start ){
  register int i;
  register site *s;
  register int dir,otherparity=0;
  
  /* allocate temporary work space only if not already allocated */
  if(temp_not_allocated)
    {
      for( dir=XUP; dir<=TUP; dir++ ){
	temp[dir]  =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      }
      temp_not_allocated = 0 ;
    }
  
  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Start gathers from positive directions */
  for(dir=XUP; dir<=TUP; dir++){
    if(start==1) tag[dir] = start_gather_field( src, sizeof(su3_vector),
					       dir, parity, gen_pt[dir] );
    else restart_gather_field( src, sizeof(su3_vector),
			       dir, parity, gen_pt[dir] , tag[dir] );
  }
  
  /* Multiply by adjoint matrix at other sites */
  FORSOMEPARITYDOMAIN(i,s,otherparity){
    if( i < loopend-FETCH_UP ){
       prefetch_V(&(src[i+FETCH_UP]));
       prefetch_4MVVVV( 
		       &((s+FETCH_UP)->link[XUP]),
		       &(temp[0][i+FETCH_UP]),
		       &(temp[1][i+FETCH_UP]),
		       &(temp[2][i+FETCH_UP]),
		       &(temp[3][i+FETCH_UP]) );
    }
    mult_adj_su3_mat_4vec( s->link, &(src[i]), &(temp[0][i]),
			   &(temp[1][i]), &(temp[2][i]), &(temp[3][i]) );
  } END_LOOP
      
  /* Start gathers from negative directions */
  for( dir=XUP; dir <= TUP; dir++){
    if (start==1) tag[OPP_DIR(dir)] = start_gather_field( temp[dir],
	  sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
    else restart_gather_field( temp[dir], sizeof(su3_vector), 
	   OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)] );
  }
  
  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[dir]);
  }
  
#ifdef SCHROED_FUN
  FORSOMEPARITY(i,s,parity) if(s->t > 0){
    if(s->t == (nt-1)){
      mult_su3_mat_vec( &(s->link[XUP]),
			(su3_vector *)(gen_pt[XUP][i]), &(dest[i]));
      for(dir=YUP; dir<TUP; dir++){
	mult_su3_mat_vec_sum( &(s->link[dir]),
			      (su3_vector *)(gen_pt[dir][i]), &(dest[i]));
      }
    }
    else{
#else
      FORSOMEPARITY(i,s,parity){
	if( i < loopend-FETCH_UP ){
	  prefetch_4MVVVV(&((s+FETCH_UP)->link[XUP]), 
			  (su3_vector *)gen_pt[XUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[YUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
			  (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
	}
#endif
	mult_su3_mat_vec_sum_4dir( s->link,
		   (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
		   (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
		   &(dest[i]));
	/*------------------------------------------------------------*/
#ifdef SCHROED_FUN
      }
#endif
  } END_LOOP
    
  /*------------------------------------------------------------*/
  /* Wait gathers from negative directions, accumulate (negative) */
  for(dir=XUP; dir<=TUP; dir++){
    wait_gather(tag[OPP_DIR(dir)]);
  }
    
  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_VVVV( 
		    (su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
		    (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
    }
    
    sub_four_su3_vecs( &(dest[i]),
		       (su3_vector *)(gen_pt[XDOWN][i]),
		       (su3_vector *)(gen_pt[YDOWN][i]),
		       (su3_vector *)(gen_pt[ZDOWN][i]),
		       (su3_vector *)(gen_pt[TDOWN][i]) ); 
      /*------------------------------------------------------------*/
  } END_LOOP
}

