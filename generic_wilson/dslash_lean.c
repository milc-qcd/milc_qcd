/*************  dslash_lean.c *******************************/
/* MIMD version 6 */
/* Version to use less memory - do things one direction at a time */

/* 7/18/01 added dslash_special - CD */
/* 1/26/00 combined with Schroedinger functional version - UMH */

/*
  dslash(F_OFFSET(psi),F_OFFSET(mp),isign,l_parity);
  Compute SUM_dirs ( 
  ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
  )
  With U = exp(iA) the operation in the continuum limit is
  
  8 + 2*isign*gamma[dir]*(partial_dir + i A_dir)
  
  Note that isign = -1 is more conventionally used in the lattice
  community, but we use isign = +1 for our Wilson fermion operator.
  
*/

#include "generic_wilson_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

void dslash(field_offset src, field_offset dest, int isign, int parity)
{
  half_wilson_vector hwv;
  
  register int i,end;
  register site *s;
  register int dir,otherparity;
  msg_tag *tag[2];
  
  switch(parity) {
  case EVEN:      otherparity=ODD; break;
  case ODD:       otherparity=EVEN; break;
  case EVENANDODD:        otherparity=EVENANDODD; break;
  }
  
#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w */
    if(MAXHTMP < 2){
      printf("dslash: MAXHTMP must be 2 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 2){
      printf("dslash: N_POINTERS must be 2 or more!\n");
      terminate(1);
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
	
	tag[1]=start_gather( F_OFFSET(htmp[1]), 
			     sizeof(half_wilson_vector), OPP_DIR(dir),
			     parity, gen_pt[1] );
    
    
    /* Take Wilson projection for src displaced in up direction, gathered,
       multiply it by link matrix, expand it, and add.
       to dest */
    wait_gather(tag[0]);
    if(dir==XUP){
      FORSOMEPARITYDOMAIN(i,s,parity){
	if( i < loopend-FETCH_UP ){
	  prefetch_M( &((s+FETCH_UP)->link[dir]) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[0][i+FETCH_UP]) );
	}
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
    else if(dir==TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0 && s->t != (nt-1)){
	if( i < loopend-FETCH_UP ){
	  prefetch_M( &((s+FETCH_UP)->link[dir]) );
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[0][i+FETCH_UP]) );
	}
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#endif
    }
    else{
      FORSOMEPARITYDOMAIN(i,s,parity){
	if( i < loopend-FETCH_UP ){
	  prefetch_M( &((s+FETCH_UP)->link[dir]) );
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[0][i+FETCH_UP]) );
	}
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
	  }
    
    cleanup_gather(tag[0]);
    
    /* Take Wilson projection for src displaced in down direction,
       expand it, and add to dest */
    wait_gather(tag[1]);
    
#ifdef SCHROED_FUN
    if(dir < TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0){
	if( i < loopend-FETCH_UP ){
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[1][i+FETCH_UP]) );
	}
	wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
    }
    else{
      FORSOMEPARITY(i,s,parity) if(s->t > 1){
#else
      FORSOMEPARITY(i,s,parity){
#endif
	if( i < loopend-FETCH_UP ){
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[1][i+FETCH_UP]) );
	}
	wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
#endif
    cleanup_gather(tag[1]);
    
  } /* end loop over directions */
} /* end (of dslash) */


/*********************************************************************/
/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=0,use
  start_gather, otherwise use restart_gather.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
/*********************************************************************/

void dslash_special(field_offset src,field_offset dest,
  int isign,int parity,msg_tag **tag,int is_started)
{
  half_wilson_vector hwv;
  
  register int i,end;
  register site *s;
  register int dir,otherparity;
  
  switch(parity) {
  case EVEN:      otherparity=ODD; break;
  case ODD:       otherparity=EVEN; break;
  case EVENANDODD:        otherparity=EVENANDODD; break;
  }
  
#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w */
    if(MAXHTMP < 2){
      printf("dslash_special: MAXHTMP must be 2 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 8){
      printf("dslash_special: N_POINTERS must be 8 or more!\n");
      terminate(1);
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
	if(is_started==0)tag[dir]=start_gather( F_OFFSET(htmp[0]), 
		sizeof(half_wilson_vector),
		dir, parity, gen_pt[dir] );
        else restart_gather (F_OFFSET(htmp[0]), 
		sizeof(half_wilson_vector),
		dir, parity, gen_pt[dir], tag[dir]);
    
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
	
	if(is_started==0)
           tag[OPP_DIR(dir)] = start_gather( F_OFFSET(htmp[1]), 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
        else restart_gather(F_OFFSET(htmp[1]), 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)], tag[OPP_DIR(dir)]);    
    /* Take Wilson projection for src displaced in up direction, gathered,
       multiply it by link matrix, expand it, and add.
       to dest */
    wait_gather(tag[dir]);

    if(dir==XUP){
      FORSOMEPARITYDOMAIN(i,s,parity){
	  if( i < loopend-FETCH_UP ){
             prefetch_M( &((s+FETCH_UP)->link[dir]) );
             prefetch_H( (half_wilson_vector *)(gen_pt[dir][i+FETCH_UP]) );
        }
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[dir][i]), &hwv ); 
	wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
    else if(dir==TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0 && s->t != (nt-1)){
        if( i < loopend-FETCH_UP ){
	    prefetch_M( &((s+FETCH_UP)->link[dir]) );
            prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	    prefetch_H( (half_wilson_vector *)(gen_pt[dir][i+FETCH_UP]) );
        }
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[dir][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#endif
    }
    else{
      FORSOMEPARITYDOMAIN(i,s,parity){
        if( i < loopend-FETCH_UP ){
          prefetch_M( &((s+FETCH_UP)->link[dir]) );
          prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
          prefetch_H( (half_wilson_vector *)(gen_pt[dir][i+FETCH_UP]) );
        }
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[dir][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
	  }
    
    /* Take Wilson projection for src displaced in down direction,
       expand it, and add to dest */
    wait_gather(tag[OPP_DIR(dir)]);
    
#ifdef SCHROED_FUN
    if(dir < TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0){
        if( i < loopend-FETCH_UP ){
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i+FETCH_UP]) );
        }
	wp_grow_add( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
    }
    else{
      FORSOMEPARITY(i,s,parity) if(s->t > 1){
#else
      FORSOMEPARITY(i,s,parity){
#endif
        if( i < loopend-FETCH_UP ){
	  prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
	  prefetch_H( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i+FETCH_UP]) );
        }
	wp_grow_add( (half_wilson_vector *)(gen_pt[OPP_DIR(dir)][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
#endif
  } /* end loop over directions */
} /* dslash_lean.c */
