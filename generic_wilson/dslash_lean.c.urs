/*************  dslash_lean.c *******************************/
/* MIMD version 5 */
/* Version to use less memory - do things one direction at a time */

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

#if defined(LOOPEND) || defined(T3E)
#ifndef EVENFIRST	/* EVENFIRST is assumed defined for revised loopend */
BOMB THE COMPILE
#endif
#undef FORSOMEPARITY
#define FORSOMEPARITY(i,s,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ), s= &(lattice[i]); \
i<loopend; i++,s++)
#define END_LOOP }
#else
#define END_LOOP 	/* define it to be nothing */
#endif

#ifdef T3E
void prefetch_vector( su3_vector * );
void prefetch_wvec( wilson_vector * );
void prefetch_hwvec( half_wilson_vector * );
void prefetch_matrix( su3_matrix * );
#endif

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
  
  for(dir=XUP;dir<=TUP;dir++){
    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
#ifdef T3E
      prefetch_wvec( (wilson_vector *)(F_PT((s+1),src)) );
#endif
      wp_shrink( (wilson_vector *)F_PT(s,src), 
		 &(s->htmp[0]), dir, isign);
    } END_LOOP
	tag[0]=start_gather( F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
			     dir, parity, gen_pt[0] );
    
    /* Take Wilson projection for src displaced in down direction,
       multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
#ifdef T3E
      prefetch_wvec( (wilson_vector *)(F_PT((s+1),src)) );
      prefetch_matrix( &(s->link[dir]) );
#endif
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
#ifndef T3E
    if(dir==XUP){
#ifdef SCHROED_FUN
      FORSOMEPARITY(i,s,parity) if(s->t > 0){
#else
      FORSOMEPARITY(i,s,parity){
#endif
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
    else if(dir==TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0 && s->t != (nt-1)){
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
#endif
    }
    else{
#ifdef SCHROED_FUN
      FORSOMEPARITY(i,s,parity) if(s->t > 0){
#else
      FORSOMEPARITY(i,s,parity){
#endif
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
      } END_LOOP
	  }
#else	/* T3E */
      switch(parity){
      case EVEN:
	i=0; end = even_sites_on_node; break;
      case ODD:
	i = even_sites_on_node; end = sites_on_node; break;
      case EVENANDODD:
	i=0; end = sites_on_node; break;
      }
      end--;
      if(dir==XUP){
	for( s= &(lattice[i]) ; i<end; i++,s++ ){
	  prefetch_matrix( &(s->link[dir]) );
	  prefetch_hwvec( (half_wilson_vector *)(gen_pt[0][i+1]) );
#ifdef SCHROED_FUN
	if(s->t > 0){
#endif
	  mult_su3_mat_hwvec( &(s->link[dir]), 
			      (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	  wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
#ifdef SCHROED_FUN
	}
	}
	if(s->t > 0){
#else
	}
#endif
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
#ifdef SCHROED_FUN
      }
    }
    else if(dir==TUP){
      for( s= &(lattice[i]) ; i<end; i++,s++ ){
	prefetch_matrix( &(s->link[dir]) );
	prefetch_wvec( (wilson_vector *)(F_PT((s+1),dest)) );
	prefetch_hwvec( (half_wilson_vector *)(gen_pt[0][i+1]) );
	if(s->t > 0 && s->t != (nt-1)){
	  mult_su3_mat_hwvec( &(s->link[dir]), 
			      (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	  wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest),
		       dir, isign);
	}
      }
      if(s->t > 0 && s->t != (nt-1)){
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest),
		     dir, isign);
      }
#endif
    }
    else{
      for( s= &(lattice[i]) ; i<end; i++,s++ ){
	prefetch_matrix( &(s->link[dir]) );
	prefetch_wvec( (wilson_vector *)(F_PT((s+1),dest)) );
	prefetch_hwvec( (half_wilson_vector *)(gen_pt[0][i+1]) );
#ifdef SCHROED_FUN
	if(s->t > 0){
#endif
	  mult_su3_mat_hwvec( &(s->link[dir]), 
			      (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
	  wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
#ifdef SCHROED_FUN
	}
      }
      if(s->t > 0){
#else
      }
#endif
      mult_su3_mat_hwvec( &(s->link[dir]), 
			  (half_wilson_vector * )(gen_pt[0][i]), &hwv ); 
      wp_grow_add( &hwv, (wilson_vector *)F_PT(s,dest), dir, isign);
#ifdef SCHROED_FUN
    }
#endif
    }
#endif	/* T3E */
    
    cleanup_gather(tag[0]);
    
    /* Take Wilson projection for src displaced in down direction,
       expand it, and add to dest */
    wait_gather(tag[1]);
    
#ifndef T3E
#ifdef SCHROED_FUN
    if(dir < TUP){
      FORSOMEPARITY(i,s,parity) if(s->t > 0){
	wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
    }
    else{
      FORSOMEPARITY(i,s,parity) if(s->t > 1){
#else
      FORSOMEPARITY(i,s,parity){
#endif
	wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      } END_LOOP
#ifdef SCHROED_FUN
    }
#endif
#else	/* T3E */
      switch(parity){
      case EVEN:
	i=0; end = even_sites_on_node; break;
      case ODD:
	i = even_sites_on_node; end = sites_on_node; break;
      case EVENANDODD:
	i=0; end = sites_on_node; break;
      }
      end--;
#ifdef SCHROED_FUN
    if(dir < TUP){
      for( s= &(lattice[i]) ; i<end; i++,s++ ){
	prefetch_wvec( (wilson_vector *)(F_PT((s+1),dest)) );
	prefetch_hwvec( (half_wilson_vector *)(gen_pt[1][i+1]) );
	if(s->t > 0){
	  wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		       (wilson_vector *)F_PT(s,dest), dir, -isign);
	}
      }
      if(s->t > 0){
	wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		     (wilson_vector *)F_PT(s,dest), dir, -isign);
      }
    }
    else{
#endif
      for( s= &(lattice[i]) ; i<end; i++,s++ ){
	prefetch_wvec( (wilson_vector *)(F_PT((s+1),dest)) );
	prefetch_hwvec( (half_wilson_vector *)(gen_pt[1][i+1]) );
#ifdef SCHROED_FUN
	if(s->t > 1){
#endif
	  wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		       (wilson_vector *)F_PT(s,dest), dir, -isign);
#ifdef SCHROED_FUN
	}
      }
      if(s->t > 1){
#else
      }
#endif
      wp_grow_add( (half_wilson_vector *)(gen_pt[1][i]),
		   (wilson_vector *)F_PT(s,dest), dir, -isign);
#ifdef SCHROED_FUN
    }
    }
#endif
#endif	/* T3E */
    
    cleanup_gather(tag[1]);
    
  } /* end loop over directions */
} /* end (of dslash_w.c) */
